// All analysis runs in the visitor's browser. No recording is sent to a server.
export const DEFAULTS = Object.freeze({
  rrMin: 300, rrMax: 2000, rrDeviation: 25, rSnr: 8,
  sbpMin: 70, sbpMax: 220, dbpMin: 35, dbpMax: 140,
  ppMin: 15, ppMax: 120, sbpDrift: 8, driftBeats: 61,
  bpDelay: 180, brsSbpStep: 1, brsRrStep: 6, brsMinR: 0.85
});

const mean = a => a.reduce((s,x)=>s+x,0)/a.length;
const quantile = (a,q) => { if(!a.length)return NaN; const b=[...a].sort((x,y)=>x-y), p=(b.length-1)*q, i=Math.floor(p); return b[i]+(b[Math.min(i+1,b.length-1)]-b[i])*(p-i); };
const median = a => quantile(a,.5);
const sd = a => a.length>1 ? Math.sqrt(a.reduce((s,x)=>s+(x-mean(a))**2,0)/(a.length-1)) : NaN;
const finite = a => a.filter(Number.isFinite);
const bound = (v,lo,hi) => Math.max(lo,Math.min(hi,v));

export function inspectText(text){
  const lines=text.split(/\r?\n/).filter(x=>x.trim());
  const sample=lines.find(x=>/\d/.test(x)&&!/^\s*(Interval|Sampling|Channel)/i.test(x))||'';
  const delimiter=sample.includes('\t')?'\t':sample.includes(';')?';':sample.includes(',')?',':null;
  const split=line=>delimiter?line.split(delimiter).map(s=>s.trim()):line.trim().split(/\s+/);
  let header=[]; const first=split(lines[0]||'');
  if(first.some(x=>!Number.isFinite(Number(x)))) header=first;
  const maxCols=Math.min(16,Math.max(...lines.slice(0,30).map(l=>split(l).length),1));
  return {delimiter,headers:Array.from({length:maxCols},(_,i)=>header[i]||`Column ${i+1}`),lineCount:lines.length};
}

export function parseSignal(text,{fs=1000,valueColumn=1,timeColumn=-1,timeUnit='s',offset=0}={}){
  const meta=inspectText(text), lines=text.split(/\r?\n/);
  const split=line=>meta.delimiter?line.split(meta.delimiter):line.trim().split(/\s+/);
  const values=[], times=[]; let skipped=0;
  for(const line of lines){
    if(!line.trim())continue;
    const cols=split(line), v=Number(cols[valueColumn]);
    const t=timeColumn<0?values.length/fs:Number(cols[timeColumn])*(timeUnit==='ms'?.001:1);
    if(!Number.isFinite(v)||!Number.isFinite(t)){skipped++;continue;}
    values.push(v); times.push(t);
  }
  if(values.length<100)throw Error('Too few numeric samples. Check the selected columns and file format.');
  const t0=times[0]; for(let i=0;i<times.length;i++)times[i]=times[i]-t0+offset;
  for(let i=1;i<times.length;i++)if(times[i]<=times[i-1])throw Error('Time values must increase. Check the time column and units.');
  const dt=timeColumn<0?1/fs:median(times.slice(1,Math.min(1000,times.length)).map((t,i)=>t-times[i]));
  if(!Number.isFinite(dt)||dt<=0)throw Error('Unable to determine sampling interval.');
  return {values:Float64Array.from(values),times:Float64Array.from(times),fs:1/dt,skipped,duration:times.at(-1)-times[0]};
}

function nearestIndex(times,t){
  let lo=0,hi=times.length; while(lo<hi){const m=(lo+hi)>>1;if(times[m]<t)lo=m+1;else hi=m;}return lo;
}

export function detectRPeaks(ecg){
  const {values:x,times,fs}=ecg, n=x.length;
  if(n<fs*5)throw Error('ECG recording must contain at least five seconds.');
  const baselineRadius=Math.max(1,Math.round(fs*.08));
  const prefix=new Float64Array(n+1);for(let i=0;i<n;i++)prefix[i+1]=prefix[i]+x[i];
  const hp=new Float64Array(n), energy=new Float64Array(n), smooth=new Float64Array(n);
  for(let i=0;i<n;i++){const a=Math.max(0,i-baselineRadius),b=Math.min(n,i+baselineRadius+1);hp[i]=x[i]-(prefix[b]-prefix[a])/(b-a);}
  for(let i=1;i<n;i++){const d=hp[i]-hp[i-1];energy[i]=d*d;}
  const w=Math.max(2,Math.round(fs*.06));let sum=0;
  for(let i=0;i<n;i++){sum+=energy[i];if(i>=w)sum-=energy[i-w];smooth[i]=sum/Math.min(i+1,w);}
  const sample=[];for(let i=0;i<n;i+=Math.max(1,Math.round(fs*.025)))sample.push(smooth[i]);
  const threshold=median(sample)+.23*(quantile(sample,.98)-median(sample));
  if(!(threshold>0))throw Error('ECG amplitude is flat; no R peaks can be detected.');
  const candidates=[];let active=false,start=0;
  for(let i=0;i<n;i++){
    if(smooth[i]>threshold&&!active){active=true;start=i;}
    if(active&&(smooth[i]<=threshold||i===n-1)){
      const end=i, center=(start+end)>>1, radius=Math.round(fs*.11);
      let best=Math.max(0,center-radius),score=0;
      for(let j=Math.max(0,center-radius);j<Math.min(n,center+radius);j++)if(Math.abs(hp[j])>score){score=Math.abs(hp[j]);best=j;}
      candidates.push({index:best,time:times[best],score});active=false;
    }
  }
  candidates.sort((a,b)=>a.index-b.index);
  const peaks=[];const refractory=.3;
  for(const c of candidates){const p=peaks.at(-1);if(p&&c.time-p.time<refractory){if(c.score>p.score)peaks[peaks.length-1]=c;}else peaks.push(c);}
  return peaks;
}

export function extractBeats(ecg,bp,peaks,settings=DEFAULTS){
  const beats=[], delay=settings.bpDelay/1000;
  for(let i=0;i<peaks.length-1;i++){
    const r=peaks[i],next=peaks[i+1], rr=(next.time-r.time)*1000;
    const start=nearestIndex(bp.times,r.time+delay),end=nearestIndex(bp.times,next.time+delay);
    let sbp=-Infinity,dbp=Infinity,sbpTime=NaN,dbpTime=NaN,sum=0,count=0;
    for(let k=start;k<end;k++){const v=bp.values[k];if(v>sbp){sbp=v;sbpTime=bp.times[k];}if(v<dbp){dbp=v;dbpTime=bp.times[k];}sum+=v;count++;}
    beats.push({i,time:r.time,rr,sbp:count?sbp:NaN,dbp:count?dbp:NaN,sbpTime,dbpTime,map:count?sum/count:NaN,pp:count?sbp-dbp:NaN,score:r.score,flags:[],accepted:false});
  }
  return beats;
}

export function applyQC(beats,settings=DEFAULTS,manualReject=new Set()){
  const c={...DEFAULTS,...settings}, scores=finite(beats.map(b=>b.score));
  const scoreMed=median(scores), scoreNoise=median(scores.map(x=>Math.abs(x-scoreMed)))||scoreMed/4;
  const snrFloor=Math.max(0,scoreMed-c.rSnr*scoreNoise);
  const rrMed=median(finite(beats.map(b=>b.rr).filter(x=>x>=c.rrMin&&x<=c.rrMax)));
  const half=Math.floor(c.driftBeats/2);
  for(let i=0;i<beats.length;i++){
    const b=beats[i],f=[];
    if(manualReject.has(i))f.push('Manual rejection');
    if(b.rr<c.rrMin||b.rr>c.rrMax)f.push('RR range');
    if(Number.isFinite(rrMed)&&Math.abs(b.rr-rrMed)>rrMed*c.rrDeviation/100)f.push('RR deviation');
    if(b.score<snrFloor)f.push('Low R-peak amplitude');
    if(!Number.isFinite(b.sbp)||!Number.isFinite(b.dbp))f.push('BP missing');
    else{
      if(b.sbp<c.sbpMin||b.sbp>c.sbpMax)f.push('SBP range');
      if(b.dbp<c.dbpMin||b.dbp>c.dbpMax)f.push('DBP range');
      if(b.pp<c.ppMin||b.pp>c.ppMax)f.push('Pulse pressure range');
      const local=finite(beats.slice(Math.max(0,i-half),Math.min(beats.length,i+half+1)).map(x=>x.sbp));
      if(local.length>=5&&Math.abs(b.sbp-median(local))>c.sbpDrift)f.push('SBP local deviation');
    }
    b.flags=f;b.accepted=f.length===0;
  }
  return beats;
}

function adjacentPairs(beats,key){const a=[];for(let i=0;i<beats.length-1;i++)if(beats[i].accepted&&beats[i+1].accepted)a.push([beats[i][key],beats[i+1][key]]);return a;}
function spectralHRV(beats){
  const clean=beats.filter(b=>b.accepted);
  if(clean.length<100||clean.length/beats.length<.95)return null;
  const times=clean.map(b=>b.time), vals=clean.map(b=>b.rr),step=.25;
  const gap=times.slice(1).map((t,i)=>t-times[i]);if(Math.max(...gap)>5)return null;
  const n=Math.floor((times.at(-1)-times[0])/step);if(n<240)return null;
  const y=new Float64Array(n);let j=0;
  for(let i=0;i<n;i++){const t=times[0]+i*step;while(j<times.length-2&&times[j+1]<t)j++;const f=bound((t-times[j])/(times[j+1]-times[j]),0,1);y[i]=vals[j]*(1-f)+vals[j+1]*f;}
  const ym=mean(y), slope=(y[n-1]-y[0])/(n-1);for(let i=0;i<n;i++)y[i]=(y[i]-ym-slope*(i-(n-1)/2))*(.5-.5*Math.cos(2*Math.PI*i/(n-1)));
  let lf=0,hf=0;const fs=4,df=fs/n;
  for(let k=1;k<=Math.floor(.4/df);k++){const f=k*df;if(f<.04)continue;let re=0,im=0;for(let i=0;i<n;i++){const phase=2*Math.PI*k*i/n;re+=y[i]*Math.cos(phase);im-=y[i]*Math.sin(phase);}const power=(re*re+im*im)*2/(fs*n*.375);if(f<.15)lf+=power*df;else hf+=power*df;}
  return {lf,hf,lfHf:hf>0?lf/hf:NaN};
}
function pearson(x,y){const xm=mean(x),ym=mean(y);let num=0,xx=0,yy=0;for(let i=0;i<x.length;i++){const a=x[i]-xm,b=y[i]-ym;num+=a*b;xx+=a*a;yy+=b*b;}return num/Math.sqrt(xx*yy);}
function slope(x,y){const xm=mean(x),ym=mean(y);let num=0,den=0;for(let i=0;i<x.length;i++){num+=(x[i]-xm)*(y[i]-ym);den+=(x[i]-xm)**2;}return num/den;}

export function sequenceBRS(beats,settings=DEFAULTS){
  const c={...DEFAULTS,...settings}, sequences=[];
  // SBP at beat i is paired with the following RR interval at beat i+1.
  for(const sign of [1,-1]){
    let i=0;
    while(i<beats.length-2){
      const valid=k=>beats[k]?.accepted&&beats[k+1]?.accepted&&beats[k+2]?.accepted&&
        sign*(beats[k+1].sbp-beats[k].sbp)>=c.brsSbpStep&&
        sign*(beats[k+2].rr-beats[k+1].rr)>=c.brsRrStep;
      if(!valid(i)){i++;continue;}
      const start=i;while(valid(i))i++;
      if(i-start<2)continue; // Three or more paired SBP and following-RR values.
      const end=i,x=[],y=[];
      for(let k=start;k<=end;k++){x.push(beats[k].sbp);y.push(beats[k+1].rr);}
      const r=pearson(x,y),s=slope(x,y);
      if(Number.isFinite(r)&&r>=c.brsMinR&&s>0)sequences.push({start,end,direction:sign>0?'up':'down',r,slope:s});
      if(i===start)i++;
    }
  }
  return {sequences,count:sequences.length,up:sequences.filter(s=>s.direction==='up').length,down:sequences.filter(s=>s.direction==='down').length,meanSlope:sequences.length?mean(sequences.map(s=>s.slope)):NaN};
}

export function metrics(beats,settings=DEFAULTS){
  const clean=beats.filter(b=>b.accepted),rr=clean.map(b=>b.rr),sbp=clean.map(b=>b.sbp),dbp=clean.map(b=>b.dbp),map=clean.map(b=>b.map),pp=clean.map(b=>b.pp);
  const rrPairs=adjacentPairs(beats,'rr'),sbpPairs=adjacentPairs(beats,'sbp');
  const diffs=rrPairs.map(([a,b])=>b-a),sbpDiffs=sbpPairs.map(([a,b])=>Math.abs(b-a));
  const hrv={meanNN:rr.length?mean(rr):NaN,meanHR:rr.length?60000/mean(rr):NaN,sdnn:sd(rr),rmssd:diffs.length?Math.sqrt(mean(diffs.map(x=>x*x))):NaN,pnn50:diffs.length?100*diffs.filter(x=>Math.abs(x)>50).length/diffs.length:NaN,...spectralHRV(beats)};
  const bpv={sbpSD:sd(sbp),sbpCV:sbp.length?100*sd(sbp)/mean(sbp):NaN,sbpARV:sbpDiffs.length?mean(sbpDiffs):NaN,dbpSD:sd(dbp),mapSD:sd(map)};
  return {n:beats.length,accepted:clean.length,rejected:beats.length-clean.length,hrv,bp:{sbp:sbp.length?mean(sbp):NaN,dbp:dbp.length?mean(dbp):NaN,map:map.length?mean(map):NaN,pp:pp.length?mean(pp):NaN},bpv,brs:sequenceBRS(beats,settings)};
}

export function analyze(ecg,bp,windowStart,windowEnd,settings=DEFAULTS,manualReject=new Set(),editedPeaks=null){
  const peaks=(editedPeaks||detectRPeaks(ecg)).filter(p=>p.time>=windowStart&&p.time<=windowEnd).sort((a,b)=>a.time-b.time);
  if(peaks.length<5)throw Error('Too few R peaks in the selected window.');
  const beats=applyQC(extractBeats(ecg,bp,peaks,settings),settings,manualReject);
  return {peaks,beats,summary:metrics(beats,settings),windowStart,windowEnd};
}

export function synthetic(duration=330,fs=250){
  const n=Math.floor(duration*fs),et=new Float64Array(n),ev=new Float64Array(n),bt=new Float64Array(n),bv=new Float64Array(n),peakTimes=[];
  let t=.5,i=0;while(t<duration-.5){peakTimes.push(t);i++;t+=.83+.045*Math.sin(i*.31)+.025*Math.sin(i*.11);}
  let p=0;
  for(let j=0;j<n;j++){
    const time=j/fs;et[j]=bt[j]=time;
    while(p<peakTimes.length-1&&peakTimes[p+1]<time)p++;
    const phase=time-peakTimes[p];
    ev[j]=.015*Math.sin(2*Math.PI*1.2*time)+1.1*Math.exp(-.5*(phase/.014)**2)-.11*Math.exp(-.5*((phase-.04)/.03)**2);
    const sbp=121+4*Math.sin(2*Math.PI*.1*peakTimes[p])+2*Math.sin(p*.19);
    bv[j]=78+(sbp-78)*Math.exp(-.5*((phase-.33)/.10)**2)+1.0*Math.sin(2*Math.PI*.03*time);
  }
  return {ecg:{times:et,values:ev,fs,skipped:0,duration},bp:{times:bt,values:bv,fs,skipped:0,duration},truth:peakTimes};
}
