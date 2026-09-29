import test from 'node:test';
import assert from 'node:assert/strict';
import {DEFAULTS,parseSignal,detectRPeaks,detectBPPulses,analyze,synthetic,sequenceBRS,metrics,applyQC,extractBeats} from '../engine.js';

test('parser skips LabChart metadata and preserves sample timing',()=>{
  const x=parseSignal('Interval= 0.01 s\n0.0,1\n0.01,2\n0.02,3\n'+Array.from({length:120},(_,i)=>`${(i+3)*.01},${i}`).join('\n'),{valueColumn:1,timeColumn:0});
  assert.equal(x.values.length,123);assert.equal(x.skipped,1);assert.ok(Math.abs(x.fs-100)<.01);assert.ok(Math.abs(x.times[1]-.01)<1e-8);
});
test('synthetic ECG peaks and BP metrics are recovered over a five-minute window',()=>{
  const s=synthetic(330,250),peaks=detectRPeaks(s.ecg),truth=s.truth.filter(t=>t>=30&&t<=330);
  const matches=truth.filter(t=>peaks.some(p=>Math.abs(p.time-t)<.04));
  assert.ok(matches.length/truth.length>.99);
  const a=analyze(s.ecg,s.bp,30,330);
  assert.ok(a.summary.accepted/a.summary.n>.98);
  assert.ok(a.summary.bp.sbp>115&&a.summary.bp.sbp<130);
  assert.ok(a.summary.bp.dbp>70&&a.summary.bp.dbp<85);
  assert.ok(a.summary.hrv.rmssd>0);
  assert.ok(Number.isFinite(a.summary.hrv.lf));
  assert.ok(a.summary.brs.count>0);
  assert.ok(a.summary.matched/a.summary.n>.98);
  assert.ok(Number.isFinite(a.summary.bpv.sbpLF));
  assert.equal(new Set(a.beats.map(b=>b.footTime)).size,a.summary.matched);
});
test('QRS detector stays on R waves with broad T waves, baseline drift, and weaker or inverted beats',()=>{
  for(const fs of [250,1000])for(const polarity of [1,-1]){
    const duration=20,n=duration*fs,times=Float64Array.from({length:n},(_,i)=>i/fs),values=new Float64Array(n),truth=[];
    for(let t=.6,k=0;t<duration-.5;t+=.78+.04*Math.sin(k++*.2))truth.push(t);
    for(let i=0;i<n;i++){
      const t=times[i];let y=.1*Math.sin(2*Math.PI*.2*t)+.03*Math.sin(2*Math.PI*27*t);
      for(let k=0;k<truth.length;k++){
        const d=t-truth[k];if(Math.abs(d)>.5)continue;
        y+=polarity*(k%5===2?.35:1.1)*Math.exp(-.5*(d/.012)**2)
          -polarity*.16*Math.exp(-.5*((d-.04)/.023)**2)
          +polarity*.65*Math.exp(-.5*((d-.27)/.085)**2);
      }
      values[i]=y;
    }
    const peaks=detectRPeaks({times,values,fs});
    const matched=truth.filter(t=>peaks.some(p=>Math.abs(p.time-t)<.045));
    assert.equal(matched.length,truth.length,`Missed R at ${fs} Hz, polarity ${polarity}`);
    assert.equal(peaks.length,truth.length,`Extra T-wave/noise detections at ${fs} Hz, polarity ${polarity}`);
  }
});
test('QRS detector rejects premature low-amplitude interbeat false peaks',()=>{
  const fs=1000,duration=22,n=duration*fs,times=Float64Array.from({length:n},(_,i)=>i/fs),values=new Float64Array(n),truth=[];
  for(let t=.7;t<duration-.7;t+=.82)truth.push(t);
  for(let i=0;i<n;i++){
    const t=times[i];let y=.08*Math.sin(2*Math.PI*.25*t)+.025*Math.sin(2*Math.PI*31*t);
    for(let k=0;k<truth.length;k++){
      const d=t-truth[k];if(Math.abs(d)>.65)continue;
      y+=1.05*Math.exp(-.5*(d/.01)**2);
      y+=.42*Math.exp(-.5*((d-.43)/.018)**2);
      y+=.35*Math.exp(-.5*((d-.28)/.08)**2);
    }
    values[i]=y;
  }
  const peaks=detectRPeaks({times,values,fs});
  const matched=truth.filter(t=>peaks.some(p=>Math.abs(p.time-t)<.04));
  assert.equal(matched.length,truth.length);
  assert.equal(peaks.length,truth.length);
});
test('QC rejects pressure spikes and adjacent HRV differences do not bridge gaps',()=>{
  const b=Array.from({length:8},(_,i)=>({i,time:i,rr:800+i*4,sbp:120,dbp:75,map:90,pp:45,score:1,accepted:true,flags:[]}));
  b[3].sbp=250;b[3].pp=175;
  applyQC(b,{...DEFAULTS,sbpDrift:20});
  assert.ok(!b[3].accepted);assert.ok(b[3].flags.includes('SBP range'));
  const m=metrics(b);assert.equal(m.accepted,7);
  assert.ok(m.hrv.rmssd<10);
});
test('overlapping sequence BRS requires three lagged pairs and estimates known slope',()=>{
  const b=[100,102,104,106,108].map((sbp,i)=>({i,time:i,rr:800+i*10,sbp,accepted:true}));
  const r=sequenceBRS(b);assert.equal(r.count,2);assert.ok(Math.abs(r.meanSlope-5)<1e-8);
  b[2].accepted=false;assert.equal(sequenceBRS(b).count,0);
});
test('BP pulse feet are detected independently and matched once to a following ECG pulse',()=>{
  const s=synthetic(20,250),peaks=detectRPeaks(s.ecg),pulses=detectBPPulses(s.bp);
  assert.ok(pulses.length>=18);assert.ok(pulses.every(p=>p.pp>=DEFAULTS.ppMin&&p.map>p.dbp&&p.sbp>p.map));
  const beats=extractBeats(s.ecg,s.bp,peaks);
  assert.equal(new Set(beats.filter(b=>Number.isFinite(b.footTime)).map(b=>b.footTime)).size,beats.filter(b=>Number.isFinite(b.footTime)).length);
  assert.ok(beats.every(b=>!Number.isFinite(b.footTime)||(b.footTime>=b.time&&b.footTime<=b.time+.5)));
});
test('BRS-only drift rejection retains BPV and missing BP does not discard ECG-only HRV',()=>{
  const b=Array.from({length:80},(_,i)=>({i,time:i,rr:800+i%3,sbp:120,dbp:80,map:93,pp:40,footTime:i+.1,score:1}));
  b[40].sbp=136;b[40].pp=56;
  applyQC(b);assert.ok(b[40].driftRejected);assert.ok(b[40].accepted);assert.ok(!b[40].brsAccepted);
  const withDrift=metrics(b);assert.ok(withDrift.bpv.sbpSD>1);
  b[5].sbp=b[5].dbp=b[5].map=b[5].pp=NaN;applyQC(b);
  assert.ok(b[5].ecgAccepted&&!b[5].accepted);
  assert.equal(metrics(b).ecgAccepted,80);
});
test('coherence-gated alpha is estimated for a coupled LF synthetic series',()=>{
  const b=Array.from({length:340},(_,i)=>({i,time:i,rr:1000+80*Math.sin(2*Math.PI*.1*(i-1)),sbp:120+10*Math.sin(2*Math.PI*.1*i),sbpDetrended:120+10*Math.sin(2*Math.PI*.1*i),dbp:80,map:94,pp:40,footTime:i+.15,accepted:true,ecgAccepted:true,brsAccepted:true}));
  const m=metrics(b);assert.ok(m.brs.coverageLF>.3);assert.ok(m.brs.alphaLF>5&&m.brs.alphaLF<12);
});
