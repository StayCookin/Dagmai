const API = "";

const strainSelect = document.getElementById("strain-select");
const runBtn = document.getElementById("run-btn");
const strainNote = document.getElementById("strain-note");
const strainMechs = document.getElementById("strain-mechs");
const resultsStatus = document.getElementById("results-status");
const resultsTable = document.getElementById("results-table");
const resultsBody = document.getElementById("results-body");
const detailPanel = document.getElementById("detail-panel");
const detailHeadline = document.getElementById("detail-headline");
const detailRationale = document.getElementById("detail-rationale");
const detailMetrics = document.getElementById("detail-metrics");
const detailCaveats = document.getElementById("detail-caveats");
const detailClose = document.getElementById("detail-close");

let strains = [];

async function fetchJSON(url, options) {
  const res = await fetch(url, options);
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`${res.status}: ${body}`);
  }
  return res.json();
}

function renderStrainInfo(strain) {
  strainNote.textContent = strain.source_note;
  strainMechs.innerHTML = "";
  if (strain.mechanisms.length === 0) {
    const tag = document.createElement("span");
    tag.className = "tag";
    tag.textContent = "no acquired resistance";
    strainMechs.appendChild(tag);
    return;
  }
  strain.mechanisms.forEach((m) => {
    const tag = document.createElement("span");
    tag.className = "tag";
    tag.textContent = m;
    strainMechs.appendChild(tag);
  });
}

function badge(text, cls) {
  return `<span class="badge ${cls}">${text}</span>`;
}

function growthBar(fraction) {
  const pct = Math.round(fraction * 100);
  const color = fraction < 0.2 ? "var(--good)" : fraction < 0.6 ? "var(--warn)" : "var(--bad)";
  return `<div class="growth-bar-wrap"><div class="growth-bar" style="width:${pct}%;background:${color}"></div></div><span style="font-size:11px;color:var(--text-dim)">${pct}%</span>`;
}

function renderResults(results) {
  resultsBody.innerHTML = "";
  results.forEach((r) => {
    const tr = document.createElement("tr");
    tr.className = "result-row";
    tr.dataset.drugA = r.drug_a.id;
    tr.dataset.drugB = r.drug_b.id;

    const nameCell = document.createElement("td");
    nameCell.innerHTML = `<strong>${r.drug_a.display_name}</strong> + <strong>${r.drug_b.display_name}</strong>`;
    tr.appendChild(nameCell);

    const growthCell = document.createElement("td");
    if (r.bliss) {
      growthCell.innerHTML = growthBar(r.bliss.growth_fraction_combo);
    } else {
      growthCell.innerHTML = '<span style="color:var(--text-dim);font-size:12px">not FBA-simulable</span>';
    }
    tr.appendChild(growthCell);

    const blissCell = document.createElement("td");
    blissCell.innerHTML = r.bliss ? badge(r.bliss.classification, r.bliss.classification) : "&mdash;";
    tr.appendChild(blissCell);

    const confCell = document.createElement("td");
    confCell.innerHTML = badge(r.explanation.confidence.replace("_", " "), r.explanation.confidence);
    tr.appendChild(confCell);

    tr.addEventListener("click", () => showDetail(r.drug_a.id, r.drug_b.id));
    resultsBody.appendChild(tr);
  });
  resultsTable.style.display = results.length ? "table" : "none";
}

async function loadStrains() {
  strains = await fetchJSON(`${API}/api/strains`);
  strainSelect.innerHTML = "";
  strains.forEach((s) => {
    const opt = document.createElement("option");
    opt.value = s.id;
    opt.textContent = s.display_name;
    strainSelect.appendChild(opt);
  });
  renderStrainInfo(strains[0]);
}

strainSelect.addEventListener("change", () => {
  const strain = strains.find((s) => s.id === strainSelect.value);
  if (strain) renderStrainInfo(strain);
});

runBtn.addEventListener("click", async () => {
  const strainId = strainSelect.value;
  runBtn.disabled = true;
  resultsStatus.style.display = "block";
  resultsStatus.textContent = "Running FBA across the panel (first run in a fresh session builds per-drug dose-response curves, ~10-20s)...";
  resultsTable.style.display = "none";
  detailPanel.style.display = "none";
  try {
    const data = await fetchJSON(`${API}/api/simulate/rank`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ strain_id: strainId }),
    });
    resultsStatus.style.display = "none";
    renderResults(data.results);
  } catch (err) {
    resultsStatus.textContent = `Error: ${err.message}`;
  } finally {
    runBtn.disabled = false;
  }
});

async function showDetail(drugA, drugB) {
  const strainId = strainSelect.value;
  detailPanel.style.display = "block";
  detailHeadline.textContent = "Loading detail...";
  detailRationale.textContent = "";
  detailMetrics.innerHTML = "";
  detailCaveats.innerHTML = "";
  detailPanel.scrollIntoView({ behavior: "smooth", block: "start" });

  try {
    const r = await fetchJSON(
      `${API}/api/simulate/pair?strain_id=${encodeURIComponent(strainId)}&drug_a=${encodeURIComponent(drugA)}&drug_b=${encodeURIComponent(drugB)}`
    );
    detailHeadline.textContent = r.explanation.headline;
    detailRationale.textContent = r.explanation.rationale;

    const metrics = [];
    if (r.bliss) {
      metrics.push(["Growth remaining (combo)", `${Math.round(r.bliss.growth_fraction_combo * 100)}%`]);
      metrics.push(["Bliss score", r.bliss.bliss_score.toFixed(3)]);
    }
    if (r.fic) {
      metrics.push(["ΣFIC index", r.fic.sigma_fic.toFixed(2)]);
      metrics.push(["FIC class", r.fic.classification]);
    }
    if (r.resistance_a) {
      metrics.push([`${r.drug_a.display_name} potency retained`, `${Math.round(r.resistance_a.multiplier * 100)}%`]);
    }
    if (r.resistance_b) {
      metrics.push([`${r.drug_b.display_name} potency retained`, `${Math.round(r.resistance_b.multiplier * 100)}%`]);
    }
    detailMetrics.innerHTML = metrics
      .map(([label, value]) => `<div class="metric"><div class="label">${label}</div><div class="value">${value}</div></div>`)
      .join("");

    detailCaveats.innerHTML = r.explanation.caveats.map((c) => `<div class="caveat">${c}</div>`).join("");
  } catch (err) {
    detailHeadline.textContent = `Error: ${err.message}`;
  }
}

detailClose.addEventListener("click", () => {
  detailPanel.style.display = "none";
});

loadStrains();
