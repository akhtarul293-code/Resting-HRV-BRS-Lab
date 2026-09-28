# Resting HRV & BRS Lab

A browser-based research prototype for resting ECG and continuous blood pressure. It lets you inspect detected R peaks and beat-level pressure values, adjust quality-control thresholds, select a window, and export an audit-friendly CSV. Recordings are read and analyzed in your browser; this project has no upload API or database.

## Important disclaimer

**Research and educational use only. This is not a medical device and must not be used for diagnosis, treatment, or clinical decisions.** The detection and QC methods are a prototype validated only with synthetic data and unit tests. They have not been benchmarked against annotated real ECG or continuous-BP recordings. Inspect R peaks, pressure calibration and alignment, rejected beats, and analysis windows before using results in a study. Spontaneous sequence BRS describes association between SBP and the following RR interval; it does not establish reflex causality. Follow your institution's ethics and data-governance requirements for real recordings.

## Use the site

1. Choose separate ECG and BP `.csv`, `.tsv`, or numeric `.txt` files. The site accepts a single amplitude/pressure column or a time plus signal column; select the correct columns in the UI. Header and metadata lines such as `Interval=` are skipped. For files without a time column, enter the sample rate (default 1,000 Hz). Select seconds or milliseconds for explicit time columns.
2. Set a BP time offset if the files have different start times. Set the BP pulse delay to associate the pressure wave with the correct R–R interval; the default is 180 ms.
3. Choose the final five minutes or a custom interval. Review the plotted ECG and BP traces and the beat flags. Click an R peak to remove it or click a missed peak to add it. You can reject individual beats in the table.
4. Download the CSV for metric summaries, settings, and a complete beat-level QC log. **Try synthetic example** loads a sample recording without any files.

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

The time axes must represent the same recording clock after the selected BP offset. Pressure values must be in **mmHg**. Files stay on the user's device; the browser processes them in memory. The downloaded CSV is saved only when the user selects Download.

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

## Host on GitHub Pages

The repository is **ready to upload without editing code**. Extract the ZIP and upload all contents, including the hidden `.github` directory, to the root of a GitHub repository. The included workflow builds, tests and deploys the site. Relative asset paths work for both project sites (`username.github.io/repository/`) and user sites.

1. In the repository, open **Settings → Pages** and set **Build and deployment → Source** to **GitHub Actions**.
2. When you have decided on site visibility, open **Actions → Deploy HRV–BRS Lab to GitHub Pages → Run workflow**. Deployment is manual by design; uploading the files alone will not publish the page.
3. GitHub will show the URL in the workflow's `github-pages` environment and under **Settings → Pages**. No code edits, secrets or API keys are required.

**Privacy:** a private GitHub repository does **not** by itself make the Pages website private. Ordinary GitHub Pages sites are generally public; private Pages access depends on eligible enterprise settings. Keep using the existing owner-only Sites deployment if access to the website must remain private. Never commit patient recordings to the repository. The browser-local processing design does not itself restrict who may visit a GitHub Pages URL. See [GitHub's Pages publishing and visibility guidance](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).

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
