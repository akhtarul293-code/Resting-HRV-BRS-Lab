# Resting HRV & BRS Lab

A browser-based research prototype for resting ECG and continuous blood pressure. It lets you inspect detected R peaks and beat-level pressure values, adjust quality-control thresholds, select a window, and export an audit-friendly CSV. Recordings are read and analyzed in your browser; this project has no upload API or database.

## Open the application

### **[Launch Resting HRV & BRS Lab](https://akhtarul293-code.github.io/Resting-HRV-BRS-Lab/)**

No installation or account is required. Recordings are processed locally in the browser tab and are not sent to this repository.

## Important disclaimer

**Research and educational use only. This is not a medical device and must not be used for diagnosis, treatment, or clinical decisions.** The detection and QC methods are a prototype validated only with synthetic data and unit tests. They have not been benchmarked against annotated real ECG or continuous-BP recordings. Inspect R peaks, pressure calibration and alignment, rejected beats, and analysis windows before using results in a study. Spontaneous sequence BRS describes association between SBP and the following RR interval; it does not establish reflex causality. Follow your institution's ethics and data-governance requirements for real recordings.

## How to use the application

1. Open the **[live application](https://akhtarul293-code.github.io/Resting-HRV-BRS-Lab/)**.
2. Select **Try synthetic example** to explore the complete workflow without choosing a recording.
3. For your own data, choose separate ECG and continuous-BP `.csv`, `.tsv`, or numeric `.txt` files. Select the correct time and signal columns. For files without a time column, enter the sample rate (default 1,000 Hz).
4. Set a BP time offset if the files have different start times. Set the BP pulse delay to associate the pressure wave with the correct R–R interval; the default is 180 ms.
5. Choose the final five minutes or a custom interval. Review the plotted signals, detected R peaks, pressure extrema and rejected beats.
6. Click a detected R peak to remove it, click a missed R wave to add it, or reject an individual beat from the table.
7. Download the CSV containing summary metrics, settings and the complete beat-level QC log.

Example ECG file:

```csv
time_s,ecg_mV
0.000,0.02
0.001,0.03
0.002,0.01
```

Example BP file:

```csv
time_s,bp_mmHg
0.000,78.0
0.001,78.1
0.002,78.2
```

The time axes must represent the same recording clock after the selected BP offset. Pressure values must be in **mmHg**. Files stay on the user's device; the browser processes them in memory. The downloaded CSV is saved only when the user selects Download. Do not commit participant recordings to this repository.

## What is calculated

| Group | Measures |
| --- | --- |
| ECG / HRV | Mean heart rate, mean NN, SDNN, RMSSD, pNN50, LF and HF power, LF/HF |
| Beat-level BP | SBP, DBP, time-averaged MAP, pulse pressure; window means |
| BPV | SBP SD, SBP coefficient of variation, SBP average real variability, DBP SD, MAP SD |
| BRS | Spontaneous sequence-method mean slope (ms/mmHg), qualifying up and down sequence counts |

Defaults include an RR range of 300–2,000 ms, local SBP median window of 61 beats and maximum local SBP deviation of 8 mmHg. All thresholds are editable. Adjacent HRV and BPV differences do not bridge rejected beats. Frequency-domain HRV is withheld if fewer than 95% of beats pass QC, a gap exceeds five seconds, or the window is too short. BRS requires at least three paired SBP / following-RR observations with matching direction and the configured step and correlation thresholds. Missing or unqualified results display as an em dash.

## Run and test locally

Install [Node.js 22 or newer](https://nodejs.org/), then in the project folder:

```sh
npm test
npm run build
python3 -m http.server 8000 -d dist
```

Open <http://localhost:8000>. No `npm install` is required: the site and tests have no third-party JavaScript dependencies. Use a local HTTP server because ES modules are restricted under `file://` in browsers.

## For developers and contributors

The application uses plain HTML, CSS and JavaScript, with no third-party runtime dependencies. Clone or download the repository, install [Node.js 22 or newer](https://nodejs.org/), and use the local commands above. Pull requests should retain the browser-local data model, visible disclaimer and synthetic validation tests.

The included `.github/workflows/pages.yml` workflow runs the tests, builds `dist/` and deploys the static files to GitHub Pages when manually started from the Actions tab.

## Project files

| File | Purpose |
| --- | --- |
| `index.html` | Page structure, controls, and visible disclaimer |
| `style.css` | Responsive styles, with system fonts and no remote font dependency |
| `app.js` | File controls, charting, beat review, and CSV export |
| `engine.js` | Parser, R-peak detector, BP extraction, QC, HRV, BPV, and BRS |
| `tests/engine.test.js` | Synthetic tests for parsing, peaks, QC, and BRS slope |
| `.github/workflows/pages.yml` | Manual GitHub Pages deployment |
| `package.json` | Dependency-free build and test commands |

The `dist/` folder is generated by `npm run build` and uploaded automatically by the workflow. It need not be committed to GitHub.

## Limitations and methods

- ECG detection uses a moving baseline and derivative-energy threshold. It is a heuristic and may fail with ectopy, poor signal quality, reversed polarity, prominent T waves, or unusual morphology. Manual edits should be independently checked.
- BP extrema are measured inside each R-to-R interval shifted by the BP pulse delay. Poor synchronization, finger-cuff calibration drift, missing data and motion can invalidate these measures. MAP is the arithmetic mean of pressure samples in that interval, not `(SBP + 2×DBP)/3`.
- HRV spectral power uses 4 Hz interpolation, detrending and a Hann periodogram. It does not identify respiratory frequency or distinguish LF/HF as a direct sympathovagal balance.
- The default ±8 mmHg SBP local-deviation rule may flag genuine physiological excursions. Review QC and document any changes to thresholds before comparing recordings.
- Synthetic validation checks the implementation on known artificial signals; it does not provide sensitivity or specificity on real annotated datasets.

Background: [1996 ESC/NASPE HRV standards](https://www.escardio.org/static-file/Escardio/Guidelines/Scientific-Statements/guidelines-Heart-Rate-Variability-FT-1996.pdf) and [spontaneous sequence BRS method](https://pubmed.ncbi.nlm.nih.gov/8119060/).
