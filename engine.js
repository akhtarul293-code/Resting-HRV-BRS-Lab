// All analysis runs in the visitor's browser. No recording is sent to a server.
export const DEFAULTS = Object.freeze({
  rrMin: 300, rrMax: 2500, rrDeviation: 25, rSnr: 8,
  sbpMin: 70, sbpMax: 220, dbpMin: 35, dbpMax: 140,
  ppMin: 10, ppMax: 120, sbpDrift: 8, driftBeats: 61,
  bpDelay: 180, bpMatchMax: 500, brsSbpStep: .5, brsRrStep: 5, brsMinR: 0.85
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
  if(!(fs>0)||!Number.isFinite(offset))throw Error('Sampling rate and time offset must be valid.');
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

// Pulse feet are found from the positive slope-sum function of low-pass BP,
// independently of ECG; ECG is used only after pulse detection for alignment.
// The DBP search is intentionally wider than the upstroke trigger so slow
// downslopes and delayed troughs are not missed.
export function detectBPPulses(bp,settings=DEFAULTS){
  const c={...DEFAULTS,...settings},{values,times,fs}=bp,n=values.length;
  if(n<fs*3)throw Error('Continuous BP recording must contain at least three seconds.');
  const smooth=new Float64Array(n),ssf=new Float64Array(n),w=Math.max(2,Math.round(.128*fs));
  const alpha=1-Math.exp(-2*Math.PI*15/fs);let acc=0;
  smooth[0]=values[0];
  for(let i=1;i<n;i++){
    smooth[i]=smooth[i-1]+alpha*(values[i]-smooth[i-1]);
    acc+=Math.max(0,smooth[i]-smooth[i-1]);
    if(i>w)acc-=Math.max(0,smooth[i-w]-smooth[i-w-1]);
    ssf[i]=acc;
  }
  const trailing=Math.max(1,Math.round(3*fs)),deque=new Int32Array(n);let head=0,tail=0;
  const onsets=[];let above=false,last=-Infinity;
  for(let i=1;i<n;i++){
    while(tail>head&&ssf[deque[tail-1]]<=ssf[i])tail--;
    deque[tail++]=i;
    while(tail>head&&deque[head]<i-trailing)head++;
    const isAbove=ssf[i]>0&&ssf[i]>.5*ssf[deque[head]];
    if(isAbove&&!above&&i-last>=c.rrMin/1000*fs){
      const lo=Math.max(0,i-Math.round(.15*fs));let foot=lo;
      for(let k=lo+1;k<=i;k++)if(smooth[k]<smooth[foot])foot=k;
      if(!onsets.length||times[foot]-times[onsets.at(-1).foot]>=c.rrMin/1000){onsets.push({foot,trigger:i});last=i;}
      else if(smooth[foot]<smooth[onsets.at(-1).foot])onsets[onsets.length-1]={foot,trigger:i};
    }
    above=isAbove;
  }
  const intervals=[];for(let i=1;i<onsets.length;i++)intervals.push(onsets[i].foot-onsets[i-1].foot);
  const medianInterval=median(intervals.length?intervals:[Math.round(fs)]);
  const pulses=[];
  for(let i=0;i<onsets.length;i++){
    const foot=onsets[i].foot,trigger=onsets[i].trigger,end=i<onsets.length-1?onsets[i+1].foot:Math.min(n-1,foot+Math.round(medianInterval));
    const peakEnd=Math.min(end,foot+Math.round(.35*fs));
    if(peakEnd<=foot+1)continue;
    let peak=foot;for(let k=foot+1;k<peakEnd;k++)if(smooth[k]>smooth[peak])peak=k;
    let dbpIndex=foot;
    const dbpLo=Math.max(0,trigger-Math.round(.25*fs)),dbpHi=Math.min(peak,trigger);
    for(let k=dbpLo;k<=dbpHi;k++)if(smooth[k]<smooth[dbpIndex])dbpIndex=k;
    let area=0,duration=0;
    for(let k=dbpIndex;k<end;k++){const dt=times[k+1]-times[k];area+=(smooth[k]+smooth[k+1])*.5*dt;duration+=dt;}
    const sbp=smooth[peak],dbp=smooth[dbpIndex],pp=sbp-dbp;
    if(pp<c.ppMin||!(duration>0))continue;
    pulses.push({time:times[foot],matchTime:Math.max(times[foot],times[peak]-.22),footIndex:foot,dbpIndex,sbpTime:times[peak],dbpTime:times[dbpIndex],sbp,dbp,map:area/duration,pp});
  }
  return pulses;
}


export function extractBeats(ecg,bp,peaks,settings=DEFAULTS){
  const c={...DEFAULTS,...settings},beats=[],pulses=detectBPPulses(bp,c);
  let cursor=0;
  for(let i=0;i<peaks.length-1;i++){
    const r=peaks[i],next=peaks[i+1], rr=(next.time-r.time)*1000;
    while(cursor<pulses.length&&(pulses[cursor].matchTime??pulses[cursor].time)<r.time)cursor++;
    let best=-1,error=Infinity;
    for(let k=cursor;k<pulses.length&&(pulses[k].matchTime??pulses[k].time)<=Math.min(next.time,r.time+c.bpMatchMax/1000);k++){
      const delta=Math.abs((pulses[k].matchTime??pulses[k].time)-r.time-c.bpDelay/1000);
      if(delta<error){best=k;error=delta;}
    }
    const pulse=best>=0?pulses[best]:null;
    if(best>=0)cursor=best+1; // Each pulse can be used at most once.
    beats.push({i,time:r.time,rr,sbp:pulse?.sbp??NaN,dbp:pulse?.dbp??NaN,sbpTime:pulse?.sbpTime??NaN,dbpTime:pulse?.dbpTime??NaN,map:pulse?.map??NaN,pp:pulse?.pp??NaN,footTime:pulse?.time??NaN,score:r.score,flags:[],accepted:false});
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
    const b=beats[i],f=[],ecgFlags=[],bpFlags=[];
    if(manualReject.has(i))ecgFlags.push('Manual rejection');
    if(b.rr<c.rrMin||b.rr>c.rrMax)ecgFlags.push('RR range');
    if(Number.isFinite(rrMed)&&Math.abs(b.rr-rrMed)>rrMed*c.rrDeviation/100)ecgFlags.push('RR deviation');
    if(b.score<snrFloor)ecgFlags.push('Low R-peak amplitude');
    if(!Number.isFinite(b.sbp)||!Number.isFinite(b.dbp))bpFlags.push('BP unmatched');
    else{
      if(b.sbp<c.sbpMin||b.sbp>c.sbpMax)bpFlags.push('SBP range');
      if(b.dbp<c.dbpMin||b.dbp>c.dbpMax)bpFlags.push('DBP range');
      if(b.pp<c.ppMin||b.pp>c.ppMax)bpFlags.push('Pulse pressure range');
    }
    f.push(...ecgFlags,...bpFlags);
    b.ecgAccepted=ecgFlags.length===0;b.accepted=f.length===0;
    b.flags=f;b.driftRejected=false;b.brsAccepted=b.accepted;b.sbpDetrended=NaN;
  }
  // Drift exclusions affect BRS only. BPV retains the raw, otherwise valid BP.
  const good=beats.filter(b=>b.accepted);
  if(good.length>=c.driftBeats){
    const overall=mean(good.map(b=>b.sbp));
    for(let i=0;i<good.length;i++){
      const local=good.slice(Math.max(0,i-half),Math.min(good.length,i+half+1)).map(b=>b.sbp);
      const baseline=median(local),b=good[i];
      b.sbpDetrended=b.sbp-baseline+overall;
      if(local.length>=5&&Math.abs(b.sbp-baseline)>c.sbpDrift){b.driftRejected=true;b.brsAccepted=false;b.flags.push('SBP drift · BRS only');}
    }
  }else for(const b of good)b.sbpDetrended=b.sbp;
  return beats;
}

function adjacentPairs(beats,key,eligibility='accepted'){const a=[];for(let i=0;i<beats.length-1;i++)if(beats[i][eligibility]&&beats[i+1][eligibility])a.push([beats[i][key],beats[i+1][key]]);return a;}
function spectralHRV(beats){
  const clean=beats.filter(b=>b.ecgAccepted);
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

// 4-Hz interpolation and 256-sample, 50%-overlap Hann Welch spectra. Never
// interpolate across a >5-second QC gap or claim coherence from one segment.
function interpolate(beats,key){
  if(beats.length<100||beats.some((b,i)=>i&&b.time-beats[i-1].time>5))return null;
  const n=Math.floor((beats.at(-1).time-beats[0].time)*4);
  if(n<384)return null;
  const out=new Float64Array(n);let j=0;
  for(let i=0;i<n;i++){
    const t=beats[0].time+i/4;
    while(j<beats.length-2&&beats[j+1].time<t)j++;
    const a=beats[j],b=beats[j+1],p=bound((t-a.time)/(b.time-a.time),0,1);
    out[i]=a[key]*(1-p)+b[key]*p;
  }
  return out;
}
function welchPair(a,b){
  if(!a||!b||a.length!==b.length)return null;
  const size=256,step=128,bins=26,ar=new Float64Array(bins),br=new Float64Array(bins),cr=new Float64Array(bins),ci=new Float64Array(bins);
  const taper=Float64Array.from({length:size},(_,i)=>.5-.5*Math.cos(2*Math.PI*i/(size-1)));
  const norm=4*taper.reduce((v,x)=>v+x*x,0);let segments=0;
  for(let start=0;start+size<=a.length;start+=step){
    const ma=mean(a.subarray(start,start+size)),mb=mean(b.subarray(start,start+size));segments++;
    for(let k=3;k<bins;k++){
      let ax=0,ay=0,bx=0,by=0;
      for(let i=0;i<size;i++){
        const theta=2*Math.PI*k*i/size,cos=Math.cos(theta)*taper[i],sin=Math.sin(theta)*taper[i];
        ax+=(a[start+i]-ma)*cos;ay-=(a[start+i]-ma)*sin;
        bx+=(b[start+i]-mb)*cos;by-=(b[start+i]-mb)*sin;
      }
      ar[k]+=2*(ax*ax+ay*ay)/norm;br[k]+=2*(bx*bx+by*by)/norm;
      cr[k]+=2*(ax*bx+ay*by)/norm;ci[k]+=2*(ay*bx-ax*by)/norm;
    }
  }
  if(segments<2)return null;
  return {ar,br,cr,ci,segments};
}
function bandResults(s,lo,hi){
  const indices=[];for(let k=3;k<26;k++){const f=k/64;if(f>=lo&&f<hi)indices.push(k);}
  if(!indices.length)return {power:NaN,alpha:NaN,coherence:NaN,coverage:0};
  let power=0,coh=0,valid=0,aa=0,bb=0;
  for(const k of indices){const a=s.ar[k]/s.segments,b=s.br[k]/s.segments;
    const c=a>0&&b>0?((s.cr[k]**2+s.ci[k]**2)/(s.ar[k]*s.br[k])):0;
    power+=b/64;coh+=c;
    if(c>=.5){valid++;aa+=a;bb+=b;}
  }
  return {power,alpha:valid&&bb>0?Math.sqrt(aa/bb):NaN,coherence:coh/indices.length,coverage:valid/indices.length};
}
function pressureSpectra(beats){
  const bp=beats.filter(b=>b.accepted),paired=beats.slice(0,-1).flatMap((b,i)=>b.brsAccepted&&beats[i+1].brsAccepted?[{...b,nextRR:beats[i+1].rr}]:[]);
  const bpSig=interpolate(bp,'sbp'),bpSpectrum=welchPair(bpSig,bpSig);
  const rrSig=interpolate(paired,'nextRR'),sbpSig=interpolate(paired,'sbpDetrended'),cross=welchPair(rrSig,sbpSig);
  const empty={alphaLF:NaN,alphaHF:NaN,coherenceLF:NaN,coherenceHF:NaN,coverageLF:0,coverageHF:0};
  const brs=cross?(()=>{const lf=bandResults(cross,.04,.15),hf=bandResults(cross,.15,.4);return {alphaLF:lf.alpha,alphaHF:hf.alpha,coherenceLF:lf.coherence,coherenceHF:hf.coherence,coverageLF:lf.coverage,coverageHF:hf.coverage};})():empty;
  const lf=bpSpectrum?bandResults(bpSpectrum,.04,.15):null,hf=bpSpectrum?bandResults(bpSpectrum,.15,.4):null;
  return {brs,bpv:{sbpLF:lf?.power??NaN,sbpHF:hf?.power??NaN,sbpLFHF:hf?.power>0?lf.power/hf.power:NaN}};
}

export function sequenceBRS(beats,settings=DEFAULTS){
  const c={...DEFAULTS,...settings}, sequences=[];
  // SBP at beat i is paired with the following RR interval at beat i+1.
  // Overlapping three-pair windows, but never across an excluded beat.
  for(let i=0;i<=beats.length-4;i++){
    if(!beats.slice(i,i+4).every(b=>b.brsAccepted??b.accepted))continue;
    const x=beats.slice(i,i+3).map(b=>Number.isFinite(b.sbpDetrended)?b.sbpDetrended:b.sbp);
    const y=beats.slice(i+1,i+4).map(b=>b.rr);
    for(const sign of [1,-1]){
      if(![0,1].every(j=>sign*(x[j+1]-x[j])>=c.brsSbpStep&&sign*(y[j+1]-y[j])>=c.brsRrStep))continue;
      const r=pearson(x,y),s=slope(x,y);
      if(Number.isFinite(r)&&r>=c.brsMinR&&s>0)sequences.push({start:i,end:i+2,direction:sign>0?'up':'down',r,slope:s});
    }
  }
  return {sequences,count:sequences.length,up:sequences.filter(s=>s.direction==='up').length,down:sequences.filter(s=>s.direction==='down').length,meanSlope:sequences.length?mean(sequences.map(s=>s.slope)):NaN};
}

export function metrics(beats,settings=DEFAULTS){
  const clean=beats.filter(b=>b.accepted),ecgClean=beats.filter(b=>b.ecgAccepted??b.accepted),rr=ecgClean.map(b=>b.rr),sbp=clean.map(b=>b.sbp),dbp=clean.map(b=>b.dbp),map=clean.map(b=>b.map),pp=clean.map(b=>b.pp);
  const rrPairs=adjacentPairs(beats,'rr','ecgAccepted'),sbpPairs=adjacentPairs(beats,'sbp');
  const diffs=rrPairs.map(([a,b])=>b-a),sbpDiffs=sbpPairs.map(([a,b])=>Math.abs(b-a));
  const hrv={meanNN:rr.length?mean(rr):NaN,meanHR:rr.length?60000/mean(rr):NaN,sdnn:sd(rr),rmssd:diffs.length?Math.sqrt(mean(diffs.map(x=>x*x))):NaN,pnn50:diffs.length?100*diffs.filter(x=>Math.abs(x)>50).length/diffs.length:NaN,...spectralHRV(beats)};
  const cv=x=>x.length&&mean(x)?100*sd(x)/mean(x):NaN;
  const spectra=pressureSpectra(beats);
  const bpv={sbpSD:sd(sbp),sbpCV:cv(sbp),sbpARV:sbpDiffs.length?mean(sbpDiffs):NaN,dbpSD:sd(dbp),dbpCV:cv(dbp),mapSD:sd(map),mapCV:cv(map),ppSD:sd(pp),ppCV:cv(pp),...spectra.bpv};
  return {n:beats.length,accepted:clean.length,rejected:beats.length-clean.length,ecgAccepted:ecgClean.length,matched:beats.filter(b=>Number.isFinite(b.footTime)).length,driftRejected:beats.filter(b=>b.driftRejected).length,hrv,bp:{sbp:sbp.length?mean(sbp):NaN,dbp:dbp.length?mean(dbp):NaN,map:map.length?mean(map):NaN,pp:pp.length?mean(pp):NaN},bpv,brs:{...sequenceBRS(beats,settings),...spectra.brs}};
}

export function analyze(ecg,bp,windowStart,windowEnd,settings=DEFAULTS,manualReject=new Set(),editedPeaks=null){
  const peaks=(editedPeaks||detectRPeaks(ecg)).filter(p=>p.time>=windowStart&&p.time<=windowEnd).sort((a,b)=>a.time-b.time);
  if(peaks.length<5)throw Error('Too few R peaks in the selected window.');
  const beats=applyQC(extractBeats(ecg,bp,peaks,settings),settings,manualReject);
  const summary=metrics(beats,settings),flags=[];
  if(windowEnd-windowStart<300)flags.push('SHORT_RECORDING');
  if(summary.matched/summary.n<.85)flags.push('LOW_ECG_BP_MATCH_RATE');
  if(summary.matched&&summary.driftRejected/summary.matched>.2)flags.push('HIGH_DRIFT_REJECTION');
  if(!(summary.hrv.meanHR>=40&&summary.hrv.meanHR<=150))flags.push('HR_OUT_OF_RANGE');
  if(summary.brs.count<5)flags.push('TOO_FEW_BRS_SEQUENCES');
  if(summary.brs.coverageLF<.3)flags.push('LOW_LF_COHERENCE_COVERAGE');
  summary.qcFlags=flags;summary.qcStatus=flags.length?'REVIEW':'OK';
  return {peaks,beats,summary,windowStart,windowEnd};
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
