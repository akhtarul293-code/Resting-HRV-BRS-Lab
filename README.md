# Resting HRV & BRS Lab

A static, browser-local research tool for paired resting ECG and continuous arterial blood-pressure recordings. Review R peaks, pressure pulse feet, systolic/diastolic measurements, beat exclusions, and a selected analysis window; export calculations and the full QC log as CSV. No server, account, API, database, or JavaScript package installation is needed to run the site.

**Live GitHub Pages URL:** [Resting HRV & BRS Lab](https://akhtarul293-code.github.io/Resting-HRV-BRS-Lab/). The URL serves whatever version has actually been deployed from that GitHub repository. Downloading this ZIP or viewing this README does **not** update the live site; publish the files and run the manual deployment workflow to put this revision online. GitHub Pages is public, so do not use it for access-controlled hosting.

## Medical and research disclaimer

**Research/education prototype, not a medical device. Do not use for diagnosis, treatment, clinical decisions, or unsupervised patient monitoring.** Synthetic tests check code behavior, not accuracy on annotated real recordings. ECG detection, pressure calibration and alignment, artifact handling, pulse transit, and BRS validity require qualified review and validation for each study. LF/HF is not a direct measure of sympathovagal balance; spontaneous BRS association is not proof of reflex causation. Follow applicable ethics approvals and data-governance requirements. Never commit participant data to this repository.

## Use the website

1. Open the [application](https://akhtarul293-code.github.io/Resting-HRV-BRS-Lab/) or run it locally (instructions below). Choose **Try synthetic example** first to inspect an artificial 5½-minute signal.
2. For your own data, select separate ECG and continuous-BP `.csv`, `.tsv`, or numeric `.txt` files. Set the signal and time columns; if there is no time column, supply the sample rate in Hz. Use the same recording clock for both; the optional BP time offset shifts BP timestamps relative to ECG. BP samples must be in mmHg.
3. Choose the last N minutes or a custom start/end in seconds. Set the expected BP-foot delay and maximum match delay. Adjust ECG, BP, drift, and BRS thresholds as appropriate to your protocol.
4. Run the analysis. Inspect ECG and BP plots, the analysis-window overlay, beat-level values, summary QC flags, and the beat table. Click the ECG plot to remove/add an R peak; use **Reject** in the table for manual exclusion.
5. Download the CSV with settings, window, metric values, QC flags, and every beat (including separately labeled ECG, BP and BRS eligibility). A dash in the UI or an empty numeric CSV cell means the result was withheld/unavailable.

The site accepts **numeric text**, not WFDB `.hea/.dat` directly. Export synchronized signal columns to CSV/TSV first, and check the exported units and clock. Example input:

```csv
time_s,ecg_mV
0.000,0.02
0.001,0.03
0.002,0.01
```

```csv
time_s,bp_mmHg
0.000,78.0
0.001,78.1
0.002,78.2
```

These three-line examples only illustrate the schema; actual analysis needs a longer continuous record. Header and other nonnumeric lines are skipped and counted. For time-column files, timestamps must increase. The analysis requires at least 30 s, but windows shorter than 300 s receive a `SHORT_RECORDING` review flag.

## Measures and QC

| Domain | Outputs | Input and exclusions |
| --- | --- | --- |
| HRV | Mean HR, mean NN, SDNN, RMSSD, pNN50, LF/HF powers and ratio | ECG-eligible RR intervals; successive differences cannot bridge rejected ECG intervals. Spectral HRV requires ≥100 intervals, ≥95% ECG acceptance, enough time, and no gap >5 s. |
| Pressure | Mean SBP, DBP, MAP, and PP | One independently detected pressure pulse matched to one preceding R peak within the configured delay. MAP is the time-weighted foot-to-foot waveform mean. |
| BPV | SD/CV of SBP, DBP, MAP, PP; SBP average real variability; SBP LF/HF powers and ratio | Matched, otherwise accepted **raw** BP beats. ARV does not cross a rejected beat. |
| Spontaneous BRS | Mean up/down sequence slope and counts; coherence-gated LF/HF alpha, LF/HF coherence and coverage | BP+ECG-eligible beats without the BRS-only rolling-median SBP drift flag. Sequence BRS pairs SBP at beat *i* with the **following** RR interval (*i+1*). |

The default RR range is 300–2,500 ms; PP ≥10 mmHg; maximum ECG→BP-foot delay 500 ms; expected foot delay 180 ms; BRS step thresholds 0.5 mmHg and 5 ms, with Pearson *r* ≥0.85. Three-pair windows overlap, so the sequence count is a count of qualifying windows, **not independent events**. A 61-beat median and 8 mmHg deviation exclude SBP outliers **from BRS only**; those otherwise valid raw pressures remain in BPV. QC status `REVIEW` includes recording duration <300 s, match rate <85%, BRS drift-rejection rate >20%, mean HR outside 40–150 bpm, <5 BRS sequences, or LF coherence coverage <30%. Review does not suppress all numeric results: interpret each displayed estimate in the context of its QC flags.

Spectral BPV/alpha use 4 Hz interpolation and 256-sample Hann Welch segments with 50% overlap, require at least two segments, and withhold results when accepted-beat gaps exceed 5 s. Alpha includes only bins with coherence ≥0.5; an insufficiently coherent band displays an unavailable alpha. The browser uses a lightweight one-pole BP low-pass smoother, a slope-sum pulse detector, and a moving-baseline derivative-energy ECG detector. These are **not** the NeuroKit/WFDB implementations in the supplied Python reference, and exact numeric parity is not asserted. See [METHODS.md](METHODS.md) for a mapping and known limitations.

## Run and test locally

Install [Node.js 22+](https://nodejs.org/) and run from the folder containing `package.json`:

```sh
npm test
npm run build
python3 -m http.server 8000 -d dist
```

Visit <http://localhost:8000>. No `npm install` is necessary. Open through HTTP rather than `file://`, which restricts JavaScript module loading. The built `dist/` contains the standalone static site.

## Put these files on GitHub

The ZIP contains the complete repository tree, including `.github/workflows/pages.yml`. **Extract it first**: GitHub's web upload does not unpack a ZIP into a repository. The reliable way to preserve the hidden `.github` directory is to clone your existing repo locally, copy the extracted files into that clone, then commit and push using Git. Do not copy any private recordings, the extracted ZIP, or local `dist/` into the repository. On GitHub, set Pages to **GitHub Actions**, then choose **Actions → Deploy HRV–BRS Lab to GitHub Pages → Run workflow**. It is deliberately manual, because deploying Pages normally makes the website publicly accessible. Confirm you are comfortable with publication before running it. If you only want private access, run the static site locally or use an access-controlled host instead.

## Repository contents

| Path | Role |
| --- | --- |
| `index.html`, `style.css`, `app.js` | Page, plots, UI, review and export |
| `engine.js` | Signal parsing, detectors, pulse matching, QC and metrics |
| `tests/engine.test.js` | Synthetic parser, ECG/BP matching, drift, HRV/BPV and BRS tests |
| `METHODS.md` | Mapping from the supplied Python protocol to this browser implementation |
| `.github/workflows/pages.yml` | Manual, test-gated GitHub Pages deployment |
| `package.json` | Dependency-free test and build commands |

All raw recordings stay in browser memory unless you deliberately download the CSV. Hosting the **code** publicly is different from hosting **data** publicly; do not upload personal recordings to GitHub. The application loads no third-party scripts or remote fonts. The background-method links below are opened only if you click them.

Background: [1996 ESC/NASPE HRV standards](https://www.escardio.org/static-file/Escardio/Guidelines/Scientific-Statements/guidelines-Heart-Rate-Variability-FT-1996.pdf) and [spontaneous sequence BRS method](https://pubmed.ncbi.nlm.nih.gov/8119060/).
