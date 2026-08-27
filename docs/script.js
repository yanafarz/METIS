/* ===========================================================================
   Metisa plana target mining — data and behaviour layer
   Owner: Chiamin. Ghaya owns index.html, style.css and assets/.

   Contract: this file only ever touches the element IDs listed in
   docs/CONTRACT.md, and only ever writes the classes listed there.

   Every DOM lookup is optional. If an element does not exist yet, that
   feature quietly does nothing instead of throwing and killing the rest of
   the page. That means this works against a half-finished index.html.

   No build step. No CDN. No dependencies. Charts are hand-rolled SVG.

   IMPORTANT: fetch() is blocked on file:// URLs, so opening the HTML by
   double-clicking will fail with a CORS error. Serve it instead:
       python -m http.server 8000
   then visit http://localhost:8000/docs/
   =========================================================================== */

'use strict';

const CONFIG = {
  dataDir: 'data/',
  files: { proteins: 'proteins.json', candidates: 'candidates.json', qc: 'qc_summary.json' },
  expectedSchemaMajor: 1,
  /* The real proteome is ~26,500 proteins. Putting that many <tr> elements in
     the DOM makes every keystroke in the search box feel broken, so render a
     window and tell the person the truth about how many matched. */
  maxRenderedRows: 400,
  searchFields: ['protein_id', 'annotation', 'gene_symbol', 'family', 'subfamily'],
  tableColumns: [
    { field: 'protein_id',    label: 'Protein ID' },
    { field: 'annotation',    label: 'Annotation' },
    { field: 'family',        label: 'Family' },
    { field: 'target_class',  label: 'Target class', badge: true },
    { field: 'hit_identity',  label: 'Identity',     format: 'percent' },
    { field: 'hit_coverage',  label: 'Coverage',     format: 'percent' },
    { field: 'hit_evalue',    label: 'E-value',      format: 'evalue' },
    { field: 'hit_bitscore',  label: 'Bit score',    format: 'number' },
    { field: 'priority_score',label: 'Score',        format: 'score' }
  ]
};

/* Display labels. The data stays snake_case; wording on screen changes here,
   never in the JSON. */
const LABELS = {
  category: {
    detoxification: 'Detoxification',
    chitin_cuticle: 'Chitin & cuticle',
    digestion_insecticide_target: 'Digestion & insecticide targets',
    hormone_immune: 'Hormones & immune',
    other: 'Other',
    unassigned: 'Unassigned',
    contaminant: 'Contaminant'
  },
  target_class: {
    lethal_rnai_target: 'Lethal RNAi target',
    insecticide_binding_site: 'Insecticide binding site',
    resistance_breaker: 'Resistance breaker',
    development_disruptor: 'Development disruptor',
    none: 'Not a target'
  },
  contamination_status: {
    clean_insect: 'Clean insect',
    non_arthropod_eukaryote: 'Other eukaryote',
    ambiguous: 'Ambiguous',
    contaminant_bacterial: 'Bacterial contaminant',
    contaminant_fungal: 'Fungal contaminant',
    no_hit: 'No BLAST hit',
    no_qualifying_hit: 'No qualifying hit',
    below_threshold: 'No qualifying hit'
  },
  confidence: { high: 'High confidence', medium: 'Medium confidence', low: 'Low confidence' },
  annotation_source: {
    blast: 'BLAST only', interpro: 'InterPro only',
    'blast+interpro': 'BLAST + InterPro', none: 'Unannotated'
  },
  owner: { chiamin: 'Chiamin', ghaya: 'Ghaya', syaf: 'Syaf', lyana: 'Lyana' }
};

const state = {
  proteins: [],
  candidates: [],
  qc: null,
  filtered: [],
  search: '',
  filters: { category: '', target_class: '', contamination_status: '', confidence: '' },
  sort: { field: 'priority_score', direction: 'desc' },
  initialised: false,
  wired: false
};

/* ---------------------------------------------------------------- utilities */

const $  = (id) => document.getElementById(id);
const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));

const show = (el) => el && el.classList.remove('is-hidden');
const hide = (el) => el && el.classList.add('is-hidden');

function label(field, value) {
  if (value == null) return '—';
  return (LABELS[field] && LABELS[field][value]) || value;
}

const fmt = {
  percent: (v) => (v == null ? '—' : v.toFixed(1) + '%'),
  number:  (v) => (v == null ? '—' : v.toLocaleString('en-GB')),
  score:   (v) => (v == null ? '—' : v.toFixed(1)),
  text:    (v) => (v == null || v === '' ? '—' : String(v)),
  /* 0.0 is a legitimate BLAST e-value: it means below double precision. */
  evalue(v) {
    if (v == null) return '—';
    if (v === 0) return '< 1e-180';
    if (v >= 0.001) return v.toFixed(3);
    return v.toExponential(1).replace('e', 'e');
  },
  aa: (v) => (v == null ? '—' : v.toLocaleString('en-GB') + ' aa')
};

function applyFormat(kind, value) {
  return fmt[kind] ? fmt[kind](value) : fmt.text(value);
}

function badge(field, value) {
  const span = document.createElement('span');
  span.className = 'badge' + (value ? ' badge--' + value : '');
  span.textContent = label(field, value);
  return span;
}

function debounce(fn, ms = 150) {
  let timer;
  return (...args) => { clearTimeout(timer); timer = setTimeout(() => fn(...args), ms); };
}

/* Sorts nulls to the end regardless of direction, so an empty column never
   pushes real data off the top of the table. */
function compare(a, b, field, direction) {
  const va = a[field], vb = b[field];
  if (va == null && vb == null) return 0;
  if (va == null) return 1;
  if (vb == null) return -1;
  let result;
  if (typeof va === 'number' && typeof vb === 'number') result = va - vb;
  else result = String(va).localeCompare(String(vb), 'en-GB', { numeric: true });
  return direction === 'desc' ? -result : result;
}

/* ------------------------------------------------------------------ loading */

async function loadJSON(filename) {
  const url = CONFIG.dataDir + filename;
  const response = await fetch(url, { cache: 'no-store' });
  if (!response.ok) throw new Error(`${filename} — server returned ${response.status}`);
  return response.json();
}

function checkSchema(filename, doc) {
  const version = (doc.meta && doc.meta.schema_version) || '';
  const major = version.split('.')[0];
  if (major !== String(CONFIG.expectedSchemaMajor)) {
    console.warn(
      `[schema] ${filename} is v${version || 'unknown'} but this script expects ` +
      `v${CONFIG.expectedSchemaMajor}.x. Field names may have changed.`
    );
  }
  return doc;
}

async function init() {
  /* Idempotent. If this somehow runs twice — a duplicated <script> tag, a
     manual init() from the console, a dev server reload — the second call is
     a no-op. Without this, event listeners double up and each sort click
     toggles direction twice, which looks like sorting being broken. */
  if (state.initialised) {
    console.warn('[init] already initialised, ignoring repeat call');
    return;
  }
  state.initialised = true;

  const loading = $('loading');
  show(loading);

  try {
    const [proteinsDoc, candidatesDoc, qcDoc] = await Promise.all([
      loadJSON(CONFIG.files.proteins),
      loadJSON(CONFIG.files.candidates),
      loadJSON(CONFIG.files.qc)
    ]);

    checkSchema(CONFIG.files.proteins, proteinsDoc);
    checkSchema(CONFIG.files.candidates, candidatesDoc);
    checkSchema(CONFIG.files.qc, qcDoc);

    state.proteins   = proteinsDoc.proteins || [];
    state.candidates = (candidatesDoc.candidates || []).slice().sort((a, b) => a.rank - b.rank);
    state.qc         = qcDoc;

    const placeholder = [proteinsDoc, candidatesDoc, qcDoc].some(d => d.meta && d.meta.is_placeholder);
    togglePlaceholderBanner(placeholder);

    renderHeadlineStats();
    renderCandidates();
    buildTableHead();
    populateFilters();
    wireEvents();
    applyFilters();
    renderCharts();

    hide(loading);
    if (loading) loading.remove();
  } catch (error) {
    console.error(error);
    hide(loading);
    showLoadError(error);
  }

  /* Handy in the browser console: __MP.proteins, __MP.filtered, etc. */
  window.__MP = state;
}

function showLoadError(error) {
  const box = $('load-error');
  const msg = $('load-error-message');
  if (msg) {
    msg.textContent =
      `Could not load the data files. ${error.message}. ` +
      `If you opened this file directly, serve it instead: run "python -m http.server 8000" ` +
      `in the repo root and visit http://localhost:8000/docs/`;
  }
  show(box);
  if (!box) alert('Data failed to load: ' + error.message);
}

function togglePlaceholderBanner(isPlaceholder) {
  const banner = $('sample-data-banner');
  if (!banner) return;
  if (isPlaceholder) {
    if (!banner.textContent.trim()) {
      banner.textContent = 'Sample data. These are placeholder results, not our findings.';
    }
    banner.hidden = false;
    show(banner);
  } else {
    banner.hidden = true;
    hide(banner);
  }
}

/* ----------------------------------------------------------- headline stats */

function renderHeadlineStats() {
  const container = $('headline-stats');
  if (!container || !state.qc) return;

  const stats = state.qc.headline_stats || [];
  const frag = document.createDocumentFragment();

  stats.forEach(stat => {
    const item = document.createElement('div');
    item.className = 'stat';

    const value = document.createElement('span');
    value.className = 'stat-value';
    value.textContent = typeof stat.value === 'number'
      ? stat.value.toLocaleString('en-GB')
      : stat.value;

    const text = document.createElement('span');
    text.className = 'stat-label';
    text.textContent = stat.label;

    if (stat.status && stat.status !== 'verified') item.dataset.status = stat.status;

    item.append(value, text);
    frag.appendChild(item);
  });

  container.replaceChildren(frag);
}

/* --------------------------------------------------------- Top-10 candidates */

function renderCandidates() {
  const container = $('candidate-cards');
  if (!container) return;

  const template = $('candidate-card-template');
  const frag = document.createDocumentFragment();

  state.candidates.forEach(cand => {
    frag.appendChild(template ? fillTemplate(template, cand) : buildCandidateCard(cand));
  });

  container.replaceChildren(frag);
}

/* Preferred path: Ghaya supplies a <template> and we fill her markup, so the
   card design stays entirely hers. */
function fillTemplate(template, record) {
  const node = template.content.firstElementChild.cloneNode(true);

  if ('proteinId' in node.dataset || node.hasAttribute('data-protein-id')) {
    node.dataset.proteinId = record.protein_id;
  }

  $$('[data-field]', node).forEach(el => {
    const field = el.dataset.field;
    const value = record[field];
    if (field === 'target_class' || field === 'contamination_status' || field === 'confidence') {
      el.textContent = label(field, value);
      if (value) el.classList.add('badge--' + value);
    } else if (field === 'hit_identity' || field === 'hit_coverage') {
      el.textContent = fmt.percent(value);
    } else if (field === 'hit_evalue') {
      el.textContent = fmt.evalue(value);
    } else if (field === 'priority_score') {
      el.textContent = fmt.score(value);
    } else if (field === 'length_aa') {
      el.textContent = fmt.aa(value);
    } else {
      el.textContent = fmt.text(value);
    }
  });

  $$('[data-list]', node).forEach(el => {
    const values = record[el.dataset.list] || [];
    el.replaceChildren(...values.map(item => {
      const li = document.createElement('li');
      li.className = 'domain-chip';
      li.textContent = typeof item === 'object' ? `${item.id} ${item.name}` : String(item);
      return li;
    }));
  });

  $$('[data-action="open-detail"]', node).forEach(btn => {
    btn.addEventListener('click', () => openDetail(record.protein_id));
  });

  return node;
}

/* Fallback for when no template exists. Uses the classes from CONTRACT.md so
   Ghaya's CSS still applies. */
function buildCandidateCard(cand) {
  const card = document.createElement('article');
  card.className = 'candidate-card';
  card.dataset.proteinId = cand.protein_id;

  const rank = document.createElement('span');
  rank.className = 'candidate-rank';
  rank.textContent = '#' + cand.rank;

  const title = document.createElement('h3');
  title.textContent = cand.gene_symbol || cand.protein_id;

  const annotation = document.createElement('p');
  annotation.className = 'candidate-annotation';
  annotation.textContent = cand.annotation;

  const badges = document.createElement('div');
  badges.className = 'candidate-badges';
  badges.append(badge('target_class', cand.target_class), badge('confidence', cand.confidence));

  const score = document.createElement('p');
  score.className = 'candidate-score';
  score.textContent = `Priority score ${fmt.score(cand.priority_score)} / 100`;

  const why = document.createElement('p');
  why.className = 'candidate-why';
  why.textContent = cand.why_selected || cand.mechanism || '';

  const domains = document.createElement('ul');
  domains.className = 'domain-list';
  (cand.domains || []).forEach(d => {
    const li = document.createElement('li');
    li.className = 'domain-chip';
    li.textContent = `${d.id} ${d.name}`;
    domains.appendChild(li);
  });

  const button = document.createElement('button');
  button.type = 'button';
  button.className = 'candidate-detail-button';
  button.textContent = 'See the evidence';
  button.addEventListener('click', () => openDetail(cand.protein_id));

  card.append(rank, title, annotation, badges, score, why, domains, button);
  return card;
}

/* ----------------------------------------------------- Protein Explorer table */

function buildTableHead() {
  const table = $('candidate-table');
  if (!table) return;
  let head = table.querySelector('thead');

  /* Only build the header if Ghaya hasn't written one. */
  if (head && head.querySelector('th')) return;
  if (!head) { head = document.createElement('thead'); table.prepend(head); }

  const row = document.createElement('tr');
  CONFIG.tableColumns.forEach(col => {
    const th = document.createElement('th');
    th.scope = 'col';
    th.dataset.sort = col.field;
    th.textContent = col.label;
    row.appendChild(th);
  });
  head.replaceChildren(row);
}

function populateFilters() {
  const specs = [
    { id: 'filter-category',             field: 'category',             any: 'All categories' },
    { id: 'filter-target-class',         field: 'target_class',         any: 'All target classes' },
    { id: 'filter-contamination',        field: 'contamination_status', any: 'All screening results' },
    { id: 'filter-confidence',           field: 'confidence',           any: 'All confidence levels' }
  ];

  specs.forEach(spec => {
    const select = $(spec.id);
    if (!select) return;

    const values = [...new Set(state.proteins.map(p => p[spec.field]).filter(Boolean))];
    values.sort((a, b) => label(spec.field, a).localeCompare(label(spec.field, b)));

    const frag = document.createDocumentFragment();
    const anyOption = document.createElement('option');
    anyOption.value = '';
    anyOption.textContent = spec.any;
    frag.appendChild(anyOption);

    values.forEach(value => {
      const count = state.proteins.filter(p => p[spec.field] === value).length;
      const option = document.createElement('option');
      option.value = value;
      option.textContent = `${label(spec.field, value)} (${count})`;
      frag.appendChild(option);
    });

    select.replaceChildren(frag);
    select.dataset.field = spec.field;
  });
}

function applyFilters() {
  const term = state.search.trim().toLowerCase();

  state.filtered = state.proteins.filter(protein => {
    for (const [field, wanted] of Object.entries(state.filters)) {
      if (wanted && protein[field] !== wanted) return false;
    }
    if (!term) return true;
    return CONFIG.searchFields.some(field => {
      const value = protein[field];
      return value != null && String(value).toLowerCase().includes(term);
    });
  });

  state.filtered.sort((a, b) => compare(a, b, state.sort.field, state.sort.direction));
  renderTable();
}

function renderTable() {
  const body = $('candidate-table-body');
  const count = $('result-count');
  const empty = $('empty-state');

  const total = state.proteins.length;
  const matched = state.filtered.length;
  const cap = CONFIG.maxRenderedRows;
  const rows = matched > cap ? state.filtered.slice(0, cap) : state.filtered;

  if (count) {
    if (matched > cap) {
      count.textContent =
        `Showing the first ${cap.toLocaleString('en-GB')} of ` +
        `${matched.toLocaleString('en-GB')} matching proteins. Search or filter to narrow it down.`;
    } else if (matched === total) {
      count.textContent = `Showing all ${total.toLocaleString('en-GB')} proteins`;
    } else {
      count.textContent =
        `Showing ${matched.toLocaleString('en-GB')} of ${total.toLocaleString('en-GB')} proteins`;
    }
  }

  if (empty) (state.filtered.length === 0 ? show : hide)(empty);

  if (!body) return;
  const frag = document.createDocumentFragment();

  rows.forEach(protein => {
    const row = document.createElement('tr');
    row.dataset.proteinId = protein.protein_id;
    row.tabIndex = 0;
    if (protein.in_top10) row.classList.add('is-top10');

    CONFIG.tableColumns.forEach(col => {
      const cell = document.createElement('td');
      const value = protein[col.field];
      if (col.badge) cell.appendChild(badge(col.field, value));
      else cell.textContent = applyFormat(col.format || 'text', value);
      if (col.format && col.format !== 'text') cell.classList.add('is-numeric');
      row.appendChild(cell);
    });

    const open = () => openDetail(protein.protein_id);
    row.addEventListener('click', open);
    row.addEventListener('keydown', e => {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(); }
    });

    frag.appendChild(row);
  });

  body.replaceChildren(frag);
  markSortedColumn();
}

function markSortedColumn() {
  $$('#candidate-table th[data-sort]').forEach(th => {
    const active = th.dataset.sort === state.sort.field;
    th.classList.toggle('is-sorted', active);
    if (active) th.setAttribute('aria-sort', state.sort.direction === 'asc' ? 'ascending' : 'descending');
    else th.removeAttribute('aria-sort');
  });
}

/* ------------------------------------------------------------- detail modal */

const DETAIL_ROWS = [
  { field: 'protein_id',              label: 'Protein ID' },
  { field: 'gene_symbol',             label: 'Gene symbol' },
  { field: 'annotation',              label: 'Annotation' },
  { field: 'family',                  label: 'Family' },
  { field: 'subfamily',               label: 'Sub-family' },
  { field: 'category',                label: 'Category',        enum: true },
  { field: 'length_aa',               label: 'Length',          format: 'aa' },
  { field: 'contamination_status',    label: 'Screening',       enum: true },
  { field: 'target_class',            label: 'Target class',    enum: true },
  { field: 'confidence',              label: 'Confidence',      enum: true },
  { field: 'annotation_source',       label: 'Evidence',        enum: true },
  { field: 'priority_score',          label: 'Priority score',  format: 'score' },
  { field: 'rank',                    label: 'Top-10 rank' },
  { field: 'hit_subject_id',          label: 'Best hit' },
  { field: 'hit_subject_name',        label: 'Hit description' },
  { field: 'hit_organism',            label: 'Hit organism' },
  { field: 'hit_identity',            label: 'Identity',        format: 'percent' },
  { field: 'hit_coverage',            label: 'Coverage',        format: 'percent' },
  { field: 'hit_evalue',              label: 'E-value',         format: 'evalue' },
  { field: 'hit_bitscore',            label: 'Bit score',       format: 'number' },
  { field: 'mechanism',               label: 'Mechanism' },
  { field: 'why_selected',            label: 'Why we picked it' },
  { field: 'delivery_strategy',       label: 'Delivery' },
  { field: 'literature_support',      label: 'Published evidence' },
  { field: 'next_step',               label: 'Next step' },
  { field: 'total_blast_hits',        label: 'BLAST hits',      format: 'number' },
  { field: 'qualifying_hits',         label: 'Passed threshold',format: 'number' },
  { field: 'hit_lineage',             label: 'Hit lineage' },
  { field: 'screen_decision',         label: 'Screen decision' },
  { field: 'screen_reason',           label: 'Screen reason' },
  { field: 'owner',                   label: 'Analysed by',     enum: true },
  { field: 'split',                   label: 'Split' },
  { field: 'notes',                   label: 'Notes' }
];

function openDetail(proteinId) {
  const modal = $('target-modal');
  const body  = $('target-modal-body');
  const title = $('target-modal-title');

  const candidate = state.candidates.find(c => c.protein_id === proteinId);
  const protein   = state.proteins.find(p => p.protein_id === proteinId);
  const record    = candidate || protein;
  if (!record) { console.warn('[detail] no record for', proteinId); return; }

  if (title) {
    title.textContent = record.gene_symbol
      ? `${record.gene_symbol} — ${record.annotation}`
      : record.annotation || record.protein_id;
  }

  if (body) {
    const template = $('target-modal-template');
    if (template) {
      body.replaceChildren(fillTemplate(template, record));
    } else {
      body.replaceChildren(buildDetailBody(record));
    }
  }

  if (modal) {
    if (typeof modal.showModal === 'function' && !modal.open) modal.showModal();
    else { modal.setAttribute('open', ''); show(modal); }
  }
}

function buildDetailBody(record) {
  const frag = document.createDocumentFragment();

  if (record.in_top10) {
    const flag = document.createElement('p');
    flag.className = 'detail-top10';
    flag.textContent = `Top-10 candidate, ranked #${record.rank}`;
    frag.appendChild(flag);
  }

  const list = document.createElement('dl');
  list.className = 'detail-list';

  DETAIL_ROWS.forEach(row => {
    const value = record[row.field];
    if (value == null || value === '') return;

    const dt = document.createElement('dt');
    dt.textContent = row.label;
    const dd = document.createElement('dd');
    if (row.enum) dd.appendChild(badge(row.field, value));
    else dd.textContent = applyFormat(row.format || 'text', value);

    list.append(dt, dd);
  });
  frag.appendChild(list);

  if (record.domains && record.domains.length) {
    frag.appendChild(sectionHeading('Domains'));
    const ul = document.createElement('ul');
    ul.className = 'domain-list';
    record.domains.forEach(d => {
      const li = document.createElement('li');
      li.className = 'domain-chip';
      li.textContent = `${d.source} ${d.id} — ${d.name}`;
      ul.appendChild(li);
    });
    frag.appendChild(ul);
  }

  if (record.score_breakdown) {
    frag.appendChild(sectionHeading('How it scored'));
    frag.appendChild(scoreBars(record.score_breakdown));
  }

  if (record.specificity) {
    frag.appendChild(sectionHeading('Specificity'));
    const spec = record.specificity;
    const p = document.createElement('p');
    const status = spec.off_target_status || 'unknown';
    p.append(badge('off_target', status));
    const note = document.createElement('span');
    note.textContent = ' ' + (spec.notes || '');
    p.appendChild(note);
    frag.appendChild(p);

    if (spec.beneficial_conflicts && spec.beneficial_conflicts.length) {
      frag.appendChild(bulletList(spec.beneficial_conflicts, 'conflict-list'));
    }
  }

  if (record.risks && record.risks.length) {
    frag.appendChild(sectionHeading('Risks'));
    frag.appendChild(bulletList(record.risks, 'risk-list'));
  }

  if (record.flags && record.flags.length) {
    frag.appendChild(sectionHeading('Flags'));
    const wrap = document.createElement('p');
    record.flags.forEach(f => wrap.appendChild(badge('flag', f)));
    frag.appendChild(wrap);
  }

  return frag;
}

function sectionHeading(text) {
  const h = document.createElement('h4');
  h.className = 'detail-heading';
  h.textContent = text;
  return h;
}

function bulletList(items, className) {
  const ul = document.createElement('ul');
  ul.className = className;
  items.forEach(item => {
    const li = document.createElement('li');
    li.textContent = item;
    ul.appendChild(li);
  });
  return ul;
}

const SCORE_LABELS = {
  essentiality: 'Essentiality',
  specificity: 'Specificity',
  rnai_evidence: 'RNAi evidence',
  annotation_confidence: 'Annotation confidence',
  druggability: 'Druggability'
};

function scoreBars(breakdown) {
  const wrap = document.createElement('div');
  wrap.className = 'score-bars';

  Object.entries(breakdown).forEach(([key, value]) => {
    const row = document.createElement('div');
    row.className = 'score-row';

    const name = document.createElement('span');
    name.className = 'score-name';
    name.textContent = SCORE_LABELS[key] || key;

    const track = document.createElement('span');
    track.className = 'score-track';
    const fill = document.createElement('span');
    fill.className = 'score-fill';
    fill.style.width = Math.max(0, Math.min(100, (value / 20) * 100)) + '%';
    track.appendChild(fill);

    const num = document.createElement('span');
    num.className = 'score-number';
    num.textContent = `${value} / 20`;

    row.append(name, track, num);
    wrap.appendChild(row);
  });

  return wrap;
}

function closeDetail() {
  const modal = $('target-modal');
  if (!modal) return;
  if (typeof modal.close === 'function' && modal.open) modal.close();
  else { modal.removeAttribute('open'); hide(modal); }
}

/* -------------------------------------------------------------------- charts
   Hand-rolled SVG. No chart library, so nothing to break if a CDN is down on
   demo day. Colours come from CSS variables so Ghaya controls the palette. */

const CHART_COLOURS = [
  'var(--chart-1, #22d3ee)', 'var(--chart-2, #a78bfa)', 'var(--chart-3, #2dd4bf)',
  'var(--chart-4, #60a5fa)', 'var(--chart-5, #f472b6)', 'var(--chart-6, #fbbf24)',
  'var(--chart-7, #34d399)'
];
const TEXT_COLOUR = 'var(--chart-text, #cbd5e1)';
const MUTED_COLOUR = 'var(--chart-muted, #64748b)';

const CHART_TARGETS = [
  { id: 'chart-contamination',    key: 'contamination_breakdown' },
  { id: 'chart-sequence-quality', key: 'sequence_quality' },
  { id: 'chart-hit-rate',         key: 'hit_rate_by_split' },
  { id: 'chart-identity',         key: 'identity_distribution' },
  { id: 'chart-organisms',        key: 'top_hit_organisms' },
  { id: 'chart-thresholds',       key: 'hit_row_thresholds' },
  { id: 'chart-target-class',     key: 'candidates_by_target_class' }
];

function renderCharts() {
  if (!state.qc || !state.qc.charts) return;

  CHART_TARGETS.forEach(target => {
    const el = $(target.id);
    const spec = state.qc.charts[target.key];
    if (!el || !spec) return;

    const rows = normaliseChartData(spec);
    if (!rows.length) {
      el.replaceChildren(chartNote('No data yet.'));
      return;
    }

    const svg = spec.type === 'donut' ? donutChart(rows) : barChart(rows, spec.unit);
    const parts = [svg];
    if (spec.caption) parts.push(chartNote(spec.caption));
    if (spec.status === 'placeholder') parts.push(chartNote('Placeholder figures.', 'is-placeholder'));
    el.replaceChildren(...parts);
  });
}

/* Handles both shapes in qc_summary.json: an array of {label, value} objects,
   or parallel bins/data arrays for the histogram. */
function normaliseChartData(spec) {
  if (Array.isArray(spec.data) && spec.data.length && typeof spec.data[0] === 'object') {
    return spec.data.filter(d => d.value != null).map(d => ({ label: d.label, value: d.value }));
  }
  if (Array.isArray(spec.bins) && Array.isArray(spec.data)) {
    return spec.bins.map((bin, i) => ({ label: bin, value: spec.data[i] }))
                    .filter(d => d.value != null);
  }
  return [];
}

function chartNote(text, extraClass) {
  const p = document.createElement('p');
  p.className = 'chart-note' + (extraClass ? ' ' + extraClass : '');
  p.textContent = text;
  return p;
}

function svgEl(name, attrs = {}) {
  const el = document.createElementNS('http://www.w3.org/2000/svg', name);
  Object.entries(attrs).forEach(([k, v]) => el.setAttribute(k, v));
  return el;
}

function barChart(rows, unit = '') {
  const rowHeight = 30, gap = 8, labelWidth = 150, valueWidth = 70, padding = 8;
  const chartWidth = 320;
  const width = labelWidth + chartWidth + valueWidth;
  const height = rows.length * (rowHeight + gap) + padding * 2;
  const max = Math.max(...rows.map(r => r.value)) || 1;

  const svg = svgEl('svg', {
    viewBox: `0 0 ${width} ${height}`,
    width: '100%', height: '100%',
    preserveAspectRatio: 'xMinYMin meet',
    role: 'img',
    'aria-label': `Bar chart, ${rows.length} categories`
  });

  rows.forEach((row, i) => {
    const y = padding + i * (rowHeight + gap);
    const barWidth = Math.max(2, (row.value / max) * chartWidth);

    const text = svgEl('text', {
      x: labelWidth - 10, y: y + rowHeight / 2 + 5,
      'text-anchor': 'end', fill: TEXT_COLOUR, 'font-size': '13'
    });
    text.textContent = row.label;

    const track = svgEl('rect', {
      x: labelWidth, y, width: chartWidth, height: rowHeight,
      rx: 4, fill: MUTED_COLOUR, opacity: '0.15'
    });

    const bar = svgEl('rect', {
      x: labelWidth, y, width: barWidth, height: rowHeight,
      rx: 4, fill: CHART_COLOURS[i % CHART_COLOURS.length]
    });

    const value = svgEl('text', {
      x: labelWidth + chartWidth + 10, y: y + rowHeight / 2 + 5,
      fill: TEXT_COLOUR, 'font-size': '13', 'font-weight': '600'
    });
    value.textContent = unit === 'percent'
      ? row.value.toFixed(1) + '%'
      : row.value.toLocaleString('en-GB');

    const title = svgEl('title');
    title.textContent = `${row.label}: ${row.value.toLocaleString('en-GB')} ${unit || ''}`.trim();
    bar.appendChild(title);

    svg.append(text, track, bar, value);
  });

  return svg;
}

function donutChart(rows) {
  const size = 240, cx = size / 2, cy = size / 2, outer = 100, inner = 62;
  const total = rows.reduce((sum, r) => sum + r.value, 0) || 1;

  const width = size + 230;
  const svg = svgEl('svg', {
    viewBox: `0 0 ${width} ${size}`,
    width: '100%', height: '100%',
    preserveAspectRatio: 'xMinYMin meet',
    role: 'img',
    'aria-label': `Donut chart of ${rows.length} categories, ${total.toLocaleString('en-GB')} total`
  });

  let angle = -Math.PI / 2;

  rows.forEach((row, i) => {
    const sweep = (row.value / total) * Math.PI * 2;
    const colour = CHART_COLOURS[i % CHART_COLOURS.length];

    /* Tiny slices would be invisible, so give them a minimum sweep. */
    const drawn = Math.max(sweep, 0.012);
    const path = svgEl('path', { d: arcPath(cx, cy, outer, inner, angle, angle + drawn), fill: colour });
    const title = svgEl('title');
    const pct = ((row.value / total) * 100).toFixed(1);
    title.textContent = `${row.label}: ${row.value.toLocaleString('en-GB')} (${pct}%)`;
    path.appendChild(title);
    svg.appendChild(path);
    angle += drawn;

    /* Legend */
    const ly = 20 + i * 26;
    svg.appendChild(svgEl('rect', { x: size + 10, y: ly - 10, width: 12, height: 12, rx: 3, fill: colour }));
    const legend = svgEl('text', { x: size + 30, y: ly, fill: TEXT_COLOUR, 'font-size': '13' });
    legend.textContent = `${row.label} — ${row.value.toLocaleString('en-GB')} (${pct}%)`;
    svg.appendChild(legend);
  });

  const totalValue = svgEl('text', {
    x: cx, y: cy - 2, 'text-anchor': 'middle',
    fill: TEXT_COLOUR, 'font-size': '22', 'font-weight': '700'
  });
  totalValue.textContent = total.toLocaleString('en-GB');
  const totalLabel = svgEl('text', {
    x: cx, y: cy + 18, 'text-anchor': 'middle', fill: MUTED_COLOUR, 'font-size': '12'
  });
  totalLabel.textContent = 'proteins';
  svg.append(totalValue, totalLabel);

  return svg;
}

function arcPath(cx, cy, outer, inner, start, end) {
  const p = (radius, angle) => [cx + radius * Math.cos(angle), cy + radius * Math.sin(angle)];
  const [x1, y1] = p(outer, start);
  const [x2, y2] = p(outer, end);
  const [x3, y3] = p(inner, end);
  const [x4, y4] = p(inner, start);
  const large = end - start > Math.PI ? 1 : 0;
  return [
    `M ${x1} ${y1}`,
    `A ${outer} ${outer} 0 ${large} 1 ${x2} ${y2}`,
    `L ${x3} ${y3}`,
    `A ${inner} ${inner} 0 ${large} 0 ${x4} ${y4}`,
    'Z'
  ].join(' ');
}

/* --------------------------------------------------------------- event wiring */

function wireEvents() {
  if (state.wired) return;
  state.wired = true;

  const search = $('protein-search');
  if (search) {
    search.addEventListener('input', debounce(e => {
      state.search = e.target.value;
      applyFilters();
    }, 150));
  }

  ['filter-category', 'filter-target-class', 'filter-contamination', 'filter-confidence']
    .forEach(id => {
      const select = $(id);
      if (!select) return;
      select.addEventListener('change', e => {
        state.filters[select.dataset.field] = e.target.value;
        applyFilters();
      });
    });

  const reset = $('filter-reset');
  if (reset) {
    reset.addEventListener('click', () => {
      state.search = '';
      Object.keys(state.filters).forEach(k => { state.filters[k] = ''; });
      if (search) search.value = '';
      $$('#filter-category, #filter-target-class, #filter-contamination, #filter-confidence')
        .forEach(s => { s.value = ''; });
      applyFilters();
    });
  }

  const table = $('candidate-table');
  if (table) {
    table.addEventListener('click', e => {
      const th = e.target.closest('th[data-sort]');
      if (!th) return;
      const field = th.dataset.sort;
      if (state.sort.field === field) {
        state.sort.direction = state.sort.direction === 'asc' ? 'desc' : 'asc';
      } else {
        state.sort.field = field;
        /* Numbers read better highest-first, text reads better A-to-Z. */
        const sample = state.proteins.find(p => p[field] != null);
        state.sort.direction = typeof (sample && sample[field]) === 'number' ? 'desc' : 'asc';
      }
      applyFilters();
    });
  }

  const close = $('target-modal-close');
  if (close) close.addEventListener('click', closeDetail);

  const modal = $('target-modal');
  if (modal) {
    /* Click the backdrop to close. */
    modal.addEventListener('click', e => { if (e.target === modal) closeDetail(); });
    /* <dialog> handles Escape natively; this covers a plain div fallback. */
    document.addEventListener('keydown', e => {
      if (e.key === 'Escape' && modal.hasAttribute('open')) closeDetail();
    });
  }
}

/* ---------------------------------------------------------------------- boot */

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', init);
} else {
  init();
}
