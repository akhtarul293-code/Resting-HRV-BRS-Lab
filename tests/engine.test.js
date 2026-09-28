import test from 'node:test';
import assert from 'node:assert/strict';
import {DEFAULTS,parseSignal,detectRPeaks,analyze,synthetic,sequenceBRS,metrics,applyQC} from '../engine.js';

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
});
test('QC rejects pressure spikes and adjacent HRV differences do not bridge gaps',()=>{
  const b=Array.from({length:8},(_,i)=>({i,time:i,rr:800+i*4,sbp:120,dbp:75,map:90,pp:45,score:1,accepted:true,flags:[]}));
  b[3].sbp=250;b[3].pp=175;
  applyQC(b,{...DEFAULTS,sbpDrift:20});
  assert.ok(!b[3].accepted);assert.ok(b[3].flags.includes('SBP range'));
  const m=metrics(b);assert.equal(m.accepted,7);
  assert.ok(m.hrv.rmssd<10);
});
test('sequence BRS requires three lagged pairs and estimates known slope',()=>{
  const b=[100,102,104,106,108].map((sbp,i)=>({i,time:i,rr:800+i*10,sbp,accepted:true}));
  const r=sequenceBRS(b);assert.equal(r.count,1);assert.ok(Math.abs(r.meanSlope-5)<1e-8);
  b[2].accepted=false;assert.equal(sequenceBRS(b).count,0);
});
