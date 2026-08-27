# Front-end contract — Ghaya and Chiamin

One rule: **we never edit the same file.**

| File | Owner |
|---|---|
| `docs/index.html` | Ghaya |
| `docs/style.css` | Ghaya |
| `docs/assets/**` | Ghaya |
| `docs/script.js` | Chiamin |
| `docs/data/**` | Chiamin |
| `backend/**` | Chiamin |
| `docs/CONTRACT.md` | Both — discuss before changing |

Ghaya owns everything about how it looks. I own everything about what it does. The contract below is the seam between them: Ghaya guarantees these element IDs exist, I guarantee my JavaScript only ever touches these IDs and only ever writes the classes listed here.

Visual direction stays as agreed in the proposal: dark navy background, cyan/teal/purple accents, rounded cards, status badges, subtle hover. That's Ghaya's call, not mine.

---

## 1. Elements Ghaya provides

Structure only. Style them however you like, and nest them inside whatever wrappers the layout needs.

### Shell

| ID | Element | Purpose |
|---|---|---|
| `#sample-data-banner` | `div` | Warning strip. I show it while the data files are still placeholders and hide it when real data lands. Start it hidden with `hidden`. |
| `#load-error` | `div` | Shown if a data file fails to load. Contains `#load-error-message`. |
| `#loading` | `div` | Spinner or skeleton. I remove it once data is in. |

### Headline stats

| ID | Element | Purpose |
|---|---|---|
| `#headline-stats` | `div` or `ul` | Empty container. I inject stat items into it. |

Each item I inject looks like:

```html
<div class="stat">
  <span class="stat-value">26,490</span>
  <span class="stat-label">Proteins in the proteome</span>
</div>
```

### Top-10 candidates

| ID | Element | Purpose |
|---|---|---|
| `#candidate-cards` | `div` | Empty container. I inject ten cards. |
| `#candidate-card-template` | `<template>` | **Preferred.** If you give me a `<template>` with the markup you want, I clone and fill it instead of writing HTML in JS. Then the card design stays 100% yours. |

If you use the template, put these `data-field` attributes on the elements I should fill:

```html
<template id="candidate-card-template">
  <article class="candidate-card" data-protein-id>
    <span data-field="rank"></span>
    <h3 data-field="gene_symbol"></h3>
    <p  data-field="annotation"></p>
    <span data-field="target_class" class="badge"></span>
    <span data-field="priority_score"></span>
    <p  data-field="why_selected"></p>
    <ul data-list="domains"></ul>
    <button data-action="open-detail">See the evidence</button>
  </article>
</template>
```

I fill anything with `data-field="X"` from field `X` of the record. `data-list="domains"` gets one `<li>` per domain. `data-action="open-detail"` gets the click handler.

### Protein Explorer

| ID | Element | Purpose |
|---|---|---|
| `#protein-search` | `input type="search"` | Free-text search over ID, annotation, gene symbol, family. |
| `#filter-category` | `select` | I populate the options from the data — leave it empty. |
| `#filter-target-class` | `select` | Same. |
| `#filter-contamination` | `select` | Same. |
| `#filter-confidence` | `select` | Same. |
| `#filter-reset` | `button` | Clears search and all filters. |
| `#result-count` | `span` | I write e.g. `Showing 12 of 30 proteins`. |
| `#candidate-table` | `table` | Needs a `<thead>` with `<th data-sort="field_name">` on each sortable column, and an empty `<tbody id="candidate-table-body">`. |
| `#candidate-table-body` | `tbody` | I inject rows. |
| `#empty-state` | `div` | Shown when filters match nothing. Keep the copy actionable: "No proteins match these filters. Clear them to start again." |

Sortable columns — put `data-sort` on the `<th>` using the exact field name:
`protein_id`, `annotation`, `family`, `hit_identity`, `hit_coverage`, `hit_evalue`, `hit_bitscore`, `priority_score`.

I add `aria-sort="ascending"` / `"descending"` to the active `<th>` and a class `is-sorted`. Style those.

### Detail modal

| ID | Element | Purpose |
|---|---|---|
| `#target-modal` | `dialog` | Use a real `<dialog>` if you can — focus trapping and Escape come free. Start closed. |
| `#target-modal-body` | `div` | Empty. I inject the detail content. |
| `#target-modal-close` | `button` | Close button. |
| `#target-modal-title` | `h2` | I write the gene symbol and annotation here. |

Same trick available: give me `#target-modal-template` with `data-field` attributes and the design stays yours.

### Charts

| ID | Element | Purpose |
|---|---|---|
| `#chart-contamination` | `div` | Screening breakdown, split 3. |
| `#chart-sequence-quality` | `div` | Short sequences, X residues, internal stops. |
| `#chart-hit-rate` | `div` | Hit rate per split. |
| `#chart-identity` | `div` | Best-hit identity histogram. |
| `#chart-organisms` | `div` | Most frequent best-hit species. |

Give each one a fixed height in CSS. Charts that size themselves from a zero-height parent render as a 1-pixel line — this is the single most common way chart layouts break.

---

## 2. Classes my JavaScript writes

I only ever add and remove these. Everything else is yours.

| Class | Where | Meaning |
|---|---|---|
| `is-hidden` | any | Hidden. Please define `.is-hidden { display: none; }`. |
| `is-active` | filter chips, sort headers | Currently applied. |
| `is-sorted` | `th` | Column the table is sorted by. |
| `is-loading` | containers | Data request in flight. |
| `badge` | injected spans | Generic badge. |
| `badge--{value}` | injected spans | Value-specific badge, e.g. `badge--lethal_rnai_target`, `badge--clean_insect`, `badge--high`. |
| `stat`, `stat-value`, `stat-label` | headline stats | |
| `candidate-card` | Top-10 cards | |
| `domain-chip` | domain lists | |

Badge classes needing colours, so you can plan the palette:

```
badge--lethal_rnai_target        badge--clean_insect              badge--high
badge--insecticide_binding_site  badge--ambiguous                 badge--medium
badge--resistance_breaker        badge--contaminant_bacterial     badge--low
badge--development_disruptor     badge--contaminant_fungal
badge--none                      badge--no_hit
                                 badge--below_threshold
                                 badge--non_arthropod_eukaryote
```

Anything without a specific rule should still look acceptable from the base `.badge` styles.

---

## 3. Quality floor

Both of us hold to this, no discussion needed:

- Works down to a phone screen. A judge may open it on their own device.
- Keyboard focus is visible everywhere.
- The table is a real `<table>`, so it can be read out and copied.
- `prefers-reduced-motion` is respected.
- Nothing loads from a CDN that could be down on demo day. Vendor any chart library into `docs/assets/`.

## 4. Changing the contract

If you need an ID that isn't here, add it to this file and tell me. Don't rename an existing one without telling me first — silent renames are how the site breaks the night before submission.
