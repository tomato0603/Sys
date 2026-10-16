/* ===================================================================
 * Gibson: the interface.
 *
 * Three blocks: the reaction table (A-2), Tm (B) and the request list for
 * Claude (C), plus the result output (Markdown copy and print) for the lab
 * notebook. A-1 (single-fragment conversion) was removed on 2026-10-10 at
 * the user's request; A-2 keeps its name so the spec sections still match.
 * Every number comes from GIBSON or PB and is printed as returned; this
 * file does not round again, because rounding in two places is how a total
 * stops matching its rows (screen spec 2).
 *
 * The request and report helpers are kept pure and exported as GIBSONUI so
 * that ui-tests.js can check the text that leaves this page.
 * =================================================================== */
var GIBSONUI = (function () {
  'use strict';

  // When this page was last assembled. Update it with every change to the
  // page, so a pasted request says which version it is about.
  var BUILT_AT = '2026-10-13 09:57';

  var REQUEST_KEY = 'gibson.requests.v1';

  function pad(n) { return (n < 10 ? '0' : '') + n; }

  // Local time as YYYY-MM-DD HH:MM.
  function stamp(date) {
    return date.getFullYear() + '-' + pad(date.getMonth() + 1) + '-'
      + pad(date.getDate()) + ' ' + pad(date.getHours()) + ':'
      + pad(date.getMinutes());
  }

  // Leading/trailing blanks and blank lines go; line breaks inside stay.
  function tidy(text) {
    return String(text === undefined || text === null ? '' : text)
      .replace(/\r\n?/g, '\n').replace(/^\s+|\s+$/g, '');
  }

  function formatRequests(entries, builtAt) {
    var lines = ['【Gibson 修正リクエスト】全 ' + entries.length + ' 件',
                 'ツール: Gibson.html（作製 ' + builtAt + '）'];
    entries.forEach(function (entry, index) {
      lines.push('', '[' + (index + 1) + '] ' + entry.at, tidy(entry.text));
    });
    return lines.join('\n');
  }

  function esc(t) {
    return String(t === undefined || t === null ? '' : t)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }

  // ===================================================================
  // Result report (Markdown, and the printed page made from it)
  // ===================================================================

  /* The Markdown is the one source: the printed page is rendered from it,
   * so what is pasted into a note and what is printed cannot disagree.
   *
   * report = { reaction, tm }, each null until that block has been
   * calculated. Blocks not yet calculated are left out. */

  function cell(v) {
    var text = v === null || v === undefined || v === '' ? '—' : String(v);
    return text.replace(/\|/g, '\\|').replace(/\n/g, ' ');
  }
  function tableLines(header, rows) {
    var line = function (cells) { return '| ' + cells.map(cell).join(' | ') + ' |'; };
    // An empty header cell stays empty; '—' would read as a missing value.
    var head = '| ' + header.map(function (h) { return h === '' ? ' ' : cell(h); })
      .join(' | ') + ' |';
    return [head, '| ' + header.map(function () { return '---'; }).join(' | ') + ' |']
      .concat(rows.map(line));
  }

  function formatReport(report, builtAt, at) {
    var lines = ['# Gibson 計算結果', '',
                 '- 出力: ' + at,
                 '- ツール: Gibson.html（作製 ' + builtAt + '）'];
    var reaction = report.reaction, tm = report.tm;

    if (reaction) {
      lines.push('', '## 反応の組み立て', '');
      lines = lines.concat(tableLines(
        ['', 'Conc ng/µL', 'Length bp', 'pmol/µL', 'target pmols', 'Volume µL', '確認用 ng'],
        reaction.rows.map(function (r) {
          return [r.name, r.concentration, r.bp, r.pmolPerUl, r.pmol, r.ul, r.ng];
        }).concat([
          ['G.A. M ×' + reaction.mixFold, '', '', '', '', reaction.mixUl, ''],
          ['miliQ', '', '', '', '', reaction.waterUl, ''],
          ['Total', '', '', '', reaction.totals.pmol, reaction.totals.ul, '']])));
      if (reaction.targets) {
        lines.push('', '### 計算確認（Vector ' + reaction.vectorNg + ' ng と倍数から計算）', '');
        lines = lines.concat(tableLines(
          ['', '目標 ng', 'Volume µL', '実際 ng', '倍数（目標）', '倍数（実際）'],
          reaction.targets.map(function (t) {
            return [t.name, t.targetNg, t.ul, t.actualNg,
                    t.targetRatio === null ? '' : '×' + t.targetRatio,
                    t.actualRatio === null ? '' : '×' + t.actualRatio];
          })));
      }
      lines.push('', '### 条件設定', '');
      lines = lines.concat(tableLines(['項目', '目安', '実際', '判定'],
        reaction.checks.map(function (c) {
          return [c.label, c.target, c.value,
                  c.ok === null ? '判定なし' : c.ok ? 'OK' : 'NG'];
        })));
    }

    if (tm) {
      lines.push('', '## Tm値の計算', '');
      lines = lines.concat(tableLines(['#', '配列', '長さ mer', 'Tm ℃', 'GC %'],
        tm.rows.map(function (r, i) {
          return [i + 1, r.sequence, r.length, r.tm, r.gcPercent];
        })));
      if (tm.diff) {
        lines.push('', 'Tm 差: ' + tm.diff.gap + ' ℃（#' + tm.diff.a + ' と #'
          + tm.diff.b + (tm.diff.spread ? '、Tm が最大と最小の組' : '') + '）');
      }
    }

    lines.push('', '## 計算式と前提', '');
    if (reaction) {
      lines.push('- 650 g/mol は二本鎖 DNA 1 bp の平均分子量（慣用値）で、'
                 + 'GC 含量が偏った断片ではわずかにずれる',
                 '- pmol/µL = 濃度 × 1000 ÷ (bp × 650)、target pmols = Volume × pmol/µL、'
                 + '確認用 ng = Volume × 濃度',
                 '- G.A. Master Mix は ' + reaction.mixFold + '× なので全量 ÷ '
                 + reaction.mixFold + '、miliQ = 全量 − Master Mix − DNA の合計',
                 '- Total の pmol は丸める前の値の和',
                 '- 丸め: pmol/µL と pmols は小数第 4 位、ng は有効数字 4 桁',
                 '- Insert / Vector は確認用 ng どうしの比');
      if (reaction.targets) {
        lines.push('- 倍数から計算する方式: 目標 ng = Vector の量 × 倍数、'
                   + 'Volume = 目標 ng ÷ Conc を小数第 2 位に丸めたもの');
      }
    }
    if (tm) {
      lines.push('- Tm (℃) = (A と T の数) × 2 + (G と C の数) × 4 − 5'
                 + '（タカラバイオの簡易式）。GC % は小数第 1 位',
                 '- 簡易式の適用範囲は 25 mer 以下。超える配列の Tm は目安');
    }
    return lines.join('\n');
  }

  /* What goes on paper: the header, the headings and the tables. The
   * formula section, lists and sentences are left to the Markdown copy
   * (request of 2026-10-10: "表部分のみ"). The header block before the first
   * section is kept so the page says when it was printed. */
  function tablesOnly(markdown) {
    var out = [], inHeader = true, skipping = false;
    markdown.split('\n').forEach(function (line) {
      if (/^## /.test(line)) { inHeader = false; skipping = line === '## 計算式と前提'; }
      if (skipping) { return; }
      if (inHeader || /^#/.test(line) || /^\|/.test(line) || line === '') { out.push(line); }
    });
    return out.join('\n').replace(/\n{3,}/g, '\n\n').replace(/\s+$/, '');
  }

  // Cells are split on pipes that are not escaped, then unescaped.
  function splitRow(line) {
    return line.replace(/^\|\s?|\s?\|$/g, '').split(/\s\|\s/)
      .map(function (c) { return c.replace(/\\\|/g, '|'); });
  }

  /* Renders only what formatReport writes: #/##/### headings, "- " lists,
   * tables and plain paragraphs. It is not a general Markdown parser. */
  function reportToHtml(markdown) {
    var html = [], lines = markdown.split('\n'), i = 0;
    while (i < lines.length) {
      var line = lines[i];
      if (/^# /.test(line)) { html.push('<h1>' + esc(line.slice(2)) + '</h1>'); i++; }
      else if (/^## /.test(line)) { html.push('<h2>' + esc(line.slice(3)) + '</h2>'); i++; }
      else if (/^### /.test(line)) { html.push('<h3>' + esc(line.slice(4)) + '</h3>'); i++; }
      else if (/^- /.test(line)) {
        var items = [];
        while (i < lines.length && /^- /.test(lines[i])) {
          items.push('<li>' + esc(lines[i].slice(2)) + '</li>'); i++;
        }
        html.push('<ul>' + items.join('') + '</ul>');
      } else if (/^\|/.test(line)) {
        var header = splitRow(line), body = [];
        i += 2;   // header and the --- line
        while (i < lines.length && /^\|/.test(lines[i])) { body.push(splitRow(lines[i])); i++; }
        html.push('<table><thead><tr>' + header.map(function (c) {
          return '<th>' + esc(c) + '</th>';
        }).join('') + '</tr></thead><tbody>' + body.map(function (row) {
          return '<tr>' + row.map(function (c) { return '<td>' + esc(c) + '</td>'; }).join('')
            + '</tr>';
        }).join('') + '</tbody></table>');
      } else if (line === '') { i++; }
      else { html.push('<p>' + esc(line) + '</p>'); i++; }
    }
    return html.join('\n');
  }

  /* One line describing a saved set of conditions, for the saved list. */
  function conditionSummary(entry) {
    var rows = entry.rows || [];
    var parts = [entry.mode === 'ratio' ? '倍数から計算' : 'Volume を入力',
                 'Insert ' + Math.max(rows.length - 1, 0) + ' 本'];
    rows.forEach(function (row, i) {
      var name = row.name || GIBSON.fragmentLabel(i, rows.length);
      var bits = [];
      if (row.concentration) { bits.push(row.concentration + ' ng/µL'); }
      if (row.bp) { bits.push(row.bp + ' bp'); }
      // In the ratio mode the vector's volume is typed, the inserts' ratios.
      var ratioRow = entry.mode === 'ratio' && i > 0;
      if (ratioRow ? row.ratio : row.ul) {
        bits.push(ratioRow ? '×' + row.ratio : row.ul + ' µL');
      }
      parts.push(name + ' ' + (bits.length ? bits.join('・') : '（空欄）'));
    });
    var primers = tidy(entry.tmText).split('\n').filter(function (l) { return tidy(l); }).length;
    if (primers) { parts.push('Tm の配列 ' + primers + ' 本'); }
    return parts.join(' / ');
  }

  /* Which section the tab bar should highlight. tops are the sections'
   * top edges relative to the viewport, in page order; line is where the
   * sticky bar ends. The section in view is the last one whose top has
   * passed the line; at the bottom of the page it is the last section,
   * which may be too short ever to reach the line. */
  function sectionInView(tops, line, atBottom) {
    if (!tops.length) { return -1; }
    if (atBottom) { return tops.length - 1; }
    var index = 0;
    tops.forEach(function (top, i) { if (top <= line) { index = i; } });
    return index;
  }

  var exported = { BUILT_AT: BUILT_AT, stamp: stamp, tidy: tidy,
                   sectionInView: sectionInView,
                   conditionSummary: conditionSummary,
                   formatRequests: formatRequests,
                   formatReport: formatReport, tablesOnly: tablesOnly,
                   reportToHtml: reportToHtml };
  if (typeof document === 'undefined') { return exported; }

  // ===================================================================
  // Helpers
  // ===================================================================

  function $(id) { return document.getElementById(id); }
  function banner(kind, title, body) {
    return '<div class="banner ' + kind + '"><strong>' + esc(title) + '</strong>'
      + (body ? esc(body) : '') + '</div>';
  }
  function value(id) { return $(id).value.trim(); }

  // Engine errors are shown word for word. Anything else is a bug in this
  // page, and is labelled as one rather than dressed up as an input problem.
  function errorBanner(error) {
    var known = error && (error.name === 'GibsonError' || error.name === 'PbError');
    return banner('bad', known ? error.message
                               : '予期しないエラー: ' + (error && error.message));
  }

  function notesHtml(notes) {
    return notes.map(function (note) { return banner('warn', note); }).join('');
  }

  // ===================================================================
  // A-2. The reaction
  // ===================================================================

  var REACTION_FIELDS = ['name', 'concentration', 'bp', 'ratio', 'ul'];
  var MAX_INSERTS = 5;     // vector + 5 = 6 fragments, the top of the guide

  /* The row heading: the vector's is fixed, an insert's is a name box.
   * Left empty, the insert keeps its default name (GIBSON.fragmentLabel),
   * shown greyed in the box. */
  function rowHeading(index, count) {
    var label = GIBSON.fragmentLabel(index, count);
    if (index === 0) { return '<th scope="row">' + label + '</th>'; }
    return '<th scope="row"><input type="text" class="name" data-field="name" placeholder="'
      + esc(label) + '" aria-label="' + esc(label) + ' の名前"></th>';
  }

  function computedCell(key) { return '<td class="num calc" data-out="' + key + '"></td>'; }
  function inputCell(field, mode) {
    return '<td class="in"><input type="text" data-field="' + field
      + '" inputmode="' + mode + '"></td>';
  }

  /* Rebuilds the fragment rows for the chosen number of inserts. Values
   * already typed are carried over row by row, so changing the count does
   * not wipe the vector. Rows start empty: no default amount (spec 2). */
  function buildReactionRows() {
    var kept = readReaction();
    var count = 1 + parseInt($('insertCount').value, 10);
    var body = $('reactionRows');
    body.innerHTML = '';
    for (var i = 0; i < count; i++) {
      var tr = document.createElement('tr');
      tr.innerHTML = rowHeading(i, count)
        + inputCell('concentration', 'decimal') + inputCell('bp', 'numeric')
        + computedCell('pmolPerUl') + computedCell('pmol')
        // The vector is the reference, so it has no ratio of its own.
        + (i === 0 ? '<td class="ratio-only num">×1</td>'
                   : inputCell('ratio', 'decimal').replace('<td class="in">',
                                                        '<td class="in ratio-only">'))
        + inputCell('ul', 'decimal') + computedCell('ng');
      if (kept[i]) {
        REACTION_FIELDS.forEach(function (field) {
          var input = tr.querySelector('[data-field="' + field + '"]');
          if (input) { input.value = kept[i][field]; }
        });
      }
      body.appendChild(tr);
    }
    applyMode();
  }

  function ratioMode() {
    return document.querySelector('input[name="amountMode"]:checked').value === 'ratio';
  }

  /* In the ratio mode the inserts' volumes are calculated, so their boxes
   * are read-only and hold the result; the vector's volume is still typed.
   * The recommendation belongs to the other mode. */
  function applyMode() {
    var on = ratioMode();
    $('plan').classList.toggle('mode-ratio', on);
    Array.prototype.forEach.call(
      document.querySelectorAll('#reactionRows [data-field="ul"]'),
      function (input, i) { input.readOnly = on && i > 0; });
    clearReaction();
  }

  function readReaction() {
    return Array.prototype.map.call($('reactionRows').rows, function (tr) {
      var fragment = {};
      REACTION_FIELDS.forEach(function (field) {
        var input = tr.querySelector('[data-field="' + field + '"]');
        fragment[field] = input ? input.value.trim() : '';
      });
      return fragment;
    });
  }

  // Blue for a number that meets the guide, red for one that does not,
  // nothing where there is no rule to judge by (ok === null / undefined).
  function verdictClass(ok) { return ok === true ? 'pass' : ok === false ? 'fail' : ''; }
  function paint(el, ok) {
    el.classList.remove('pass', 'fail');
    var cls = verdictClass(ok);
    if (cls) { el.classList.add(cls); }
  }

  // Derived cells are emptied whenever the inputs no longer match them.
  function clearReaction() {
    Array.prototype.forEach.call(
      document.querySelectorAll('#reactionTable [data-out]'),
      function (td) { td.textContent = ''; paint(td, null); });
    Array.prototype.forEach.call(
      document.querySelectorAll('#reactionRows [data-field="ul"]'),
      function (input) { paint(input, null); });
    $('pResult').innerHTML = '';
    $('recResult').innerHTML = '';
    recommended = null;
    setReport('reaction', null);
  }

  function checksHtml(checks) {
    return '<div class="scroll"><table><tr><th>項目</th><th>目安</th>'
      + '<th>実際</th><th>判定</th></tr>' + checks.map(function (c) {
        var mark = c.ok === null ? '— 判定なし' : c.ok ? 'OK' : 'NG';
        var cls = verdictClass(c.ok);
        return '<tr><td>' + esc(c.label) + '</td><td>' + esc(c.target) + '</td>'
          + '<td class="num ' + cls + '">' + esc(c.value) + '</td>'
          + '<td class="' + cls + '">' + mark + '</td></tr>';
      }).join('') + '</table></div>';
  }

  function guideHtml() {
    return '<table><tr><th>項目</th><th>目安</th></tr>' + GIBSON.guide().map(function (g) {
      return '<tr><td>' + esc(g.label) + '</td><td>' + esc(g.target) + '</td></tr>';
    }).join('') + '</table>';
  }

  // The ng cell is judged by the vector amount for the vector and by the
  // ratio for an insert, since that is what the insert's ng is checked for.
  function ngOk(row, index) { return index === 0 ? row.ngOk : row.ratioOk; }

  // Target against actual ng, so the ratio-mode arithmetic can be checked.
  function targetsHtml(targets) {
    return '<h3>計算確認</h3><div class="scroll"><table><tr><th></th><th>目標 ng</th>'
      + '<th>Volume µL</th><th>実際 ng <small>Volume × Conc</small></th>'
      + '<th>倍数（目標）</th><th>倍数（実際）</th></tr>'
      + targets.map(function (t) {
          return '<tr><th scope="row">' + esc(t.name) + '</th><td class="num">' + t.targetNg
            + '</td><td class="num">' + t.ul + '</td><td class="num">' + t.actualNg
            + '</td><td class="num">' + (t.targetRatio === null ? '' : '×' + t.targetRatio)
            + '</td><td class="num">' + (t.actualRatio === null ? '' : '×' + t.actualRatio)
            + '</td></tr>';
        }).join('') + '</table></div>'
      + '<p class="hint">Volume は 0.01 µL 単位に丸めるので、実際の ng は目標から少しずれます。'
      + '判定は実際の ng で行います。</p>';
  }

  function calcReaction() {
    clearReaction();
    try {
      var result = ratioMode()
        ? GIBSON.amountsFromRatio({ fragments: readReaction() })
        : GIBSON.reaction({ fragments: readReaction() });
      Array.prototype.forEach.call($('reactionRows').rows, function (tr, i) {
        var ulInput = tr.querySelector('[data-field="ul"]');
        if (result.targets) { ulInput.value = result.rows[i].ul; }
        paint(ulInput, result.rows[i].ulOk);
        paint(tr.querySelector('[data-out="ng"]'), ngOk(result.rows[i], i));
        ['pmolPerUl', 'pmol', 'ng'].forEach(function (key) {
          tr.querySelector('[data-out="' + key + '"]').textContent = result.rows[i][key];
        });
      });
      var out = function (key, v) {
        document.querySelector('#reactionTable [data-out="' + key + '"]').textContent = v;
      };
      out('mixUl', result.mixUl);
      out('waterUl', result.waterUl);
      out('totalPmol', result.totals.pmol);
      paint(document.querySelector('#reactionTable [data-out="totalPmol"]'),
            result.totals.pmolOk);
      $('pResult').innerHTML = (result.targets ? targetsHtml(result.targets) : '')
        + checksHtml(result.checks) + notesHtml(result.notes);
      result.mixFold = GIBSON.MASTER_MIX_FOLD;
      setReport('reaction', result);
    } catch (error) {
      $('pResult').innerHTML = errorBanner(error);
    }
  }

  // The recommendation waiting to be copied into the table, or null.
  var recommended = null;

  /* Shows amounts that meet the guide. Only Conc, Length and the total are
   * needed; Amount may be empty. Nothing in the table changes until the
   * user presses the button under the recommendation. */
  function recommendAmounts() {
    clearReaction();
    var box = $('recResult');
    try {
      var answer = GIBSON.recommend({ fragments: readReaction() });
      // When the ratios had to be lowered, the colours below say which
      // numbers miss the guide; the long explanation was dropped (2026-10-10).
      var html = answer.relaxed ? ''
        : banner(answer.satisfied ? 'ok' : 'bad', answer.message);
      if (answer.result) {
        recommended = answer.result.rows.map(function (row) { return row.ul; });
        var r = answer.result;
        html += '<div class="scroll"><table><tr><th></th><th>Volume µL</th><th>ng</th>'
          + '<th>pmols</th><th>Insert / Vector</th></tr>'
          + r.rows.map(function (row, i) {
              return '<tr><th scope="row">' + esc(row.name) + '</th>'
                + '<td class="num ' + verdictClass(row.ulOk) + '">' + row.ul + '</td>'
                + '<td class="num ' + (i === 0 ? verdictClass(row.ngOk) : '') + '">'
                + row.ng + '</td><td class="num">' + row.pmol + '</td>'
                + '<td class="num ' + verdictClass(row.ratioOk) + '">'
                + (row.ratio === undefined ? '' : '×' + row.ratio) + '</td></tr>';
            }).join('')
          + '<tr class="total"><th scope="row">Total</th><td></td><td></td>'
          + '<td class="num ' + verdictClass(r.totals.pmolOk) + '">'
          + r.totals.pmol + '</td><td></td></tr></table></div>'
          + '<p class="hint">miliQ は ' + r.waterUl + ' µL になります。</p>'
          + notesHtml(r.notes)
          + '<div class="actions"><button type="button" class="primary" id="applyRec">'
          + 'この Volume を表に入れる</button></div>';
      }
      box.innerHTML = html;
      if (answer.result) { $('applyRec').addEventListener('click', applyRecommendation); }
    } catch (error) {
      box.innerHTML = errorBanner(error);
    }
  }

  function applyRecommendation() {
    if (!recommended) { return; }
    var amounts = recommended;
    Array.prototype.forEach.call($('reactionRows').rows, function (tr, i) {
      tr.querySelector('[data-field="ul"]').value = amounts[i];
    });
    calcReaction();
    $('recResult').innerHTML = banner('ok', 'おすすめの Volume を表に入れて計算しました。');
  }

  // ===================================================================
  // B. Tm
  // ===================================================================

  /* Which two rows the Tm difference is about. Two rows: those two. More:
   * the highest and the lowest Tm, so the number is the spread of the set. */
  function spreadPair(rows) {
    if (rows.length === 2) { return [0, 1]; }
    var hi = 0, lo = 0;
    rows.forEach(function (row, i) {
      if (row.tm > rows[hi].tm) { hi = i; }
      if (row.tm < rows[lo].tm) { lo = i; }
    });
    if (hi === lo) { lo = hi === 0 ? 1 : 0; }   // all equal: the gap is 0 anyway
    return hi < lo ? [hi, lo] : [lo, hi];
  }

  function calcTm() {
    var out = $('tmResult');
    out.innerHTML = '';
    setReport('tm', null);
    try {
      var rows = GIBSON.tmRows($('tmText').value);
      if (!rows.length) {
        out.innerHTML = banner('bad', '配列を入力してください。');
        return;
      }
      var html = '<div class="scroll"><table><tr><th>#</th><th>配列</th>'
        + '<th>長さ mer</th><th>Tm ℃</th><th>GC %</th><th>注意</th></tr>';
      rows.forEach(function (row, i) {
        html += '<tr' + (row.formulaApplies ? '' : ' class="outside"') + '>'
          + '<td class="num">' + (i + 1) + '</td>'
          + '<td class="seq">' + esc(row.sequence) + '</td>'
          + '<td class="num">' + row.length + '</td>'
          + '<td class="num">' + row.tm + '</td>'
          + '<td class="num">' + row.gcPercent + '</td>'
          + '<td>' + esc(row.note || '') + '</td></tr>';
      });
      html += '</table></div>';
      var diff = null;
      if (rows.length >= 2) {
        var pair = spreadPair(rows);
        var gap = GIBSON.tmDifference(rows[pair[0]].sequence, rows[pair[1]].sequence);
        diff = { gap: gap, a: pair[0] + 1, b: pair[1] + 1, spread: rows.length > 2 };
        html += banner('ok', 'Tm 差: ' + gap + ' ℃',
          '#' + diff.a + ' と #' + diff.b
          + (diff.spread ? '（Tm が最大と最小の組）' : ''));
      }
      out.innerHTML = html;
      setReport('tm', { rows: rows, diff: diff });
    } catch (error) {
      out.innerHTML = errorBanner(error);
    }
  }

  // ===================================================================
  // C. Requests
  // ===================================================================

  var requests = [];
  var storageOk = true;

  function storageWarning() {
    return banner('warn', 'この画面を閉じるとリクエストは消えます'
      + '（ブラウザに保存できませんでした）。');
  }

  /* If the saved list cannot be read, it is left in place and nothing is
   * written over it, so a glitch does not wipe requests the user kept. */
  function loadRequests() {
    try {
      var raw = window.localStorage.getItem(REQUEST_KEY);
      if (raw === null) { return; }
      var list = JSON.parse(raw);
      if (!Array.isArray(list)) { throw new Error('saved requests are not a list'); }
      requests = list.filter(function (entry) {
        return entry && typeof entry.at === 'string' && typeof entry.text === 'string';
      });
    } catch (error) {
      storageOk = false;
      if (window.console) { console.error('Gibson: requests not loaded', error); }
    }
  }

  function saveRequests() {
    if (!storageOk) { return; }
    try {
      window.localStorage.setItem(REQUEST_KEY, JSON.stringify(requests));
    } catch (error) {
      storageOk = false;
      if (window.console) { console.error('Gibson: requests not saved', error); }
    }
  }

  function renderRequests() {
    $('requestList').innerHTML = requests.length
      ? requests.map(function (entry, i) {
          return '<div class="request"><div class="when">[' + (i + 1) + '] '
            + esc(entry.at) + '</div><div class="body">' + esc(entry.text)
            + '</div></div>';
        }).join('')
      : '<p class="hint">リクエストはまだありません。</p>';
    $('copyRequests').disabled = !requests.length;
    $('clearRequests').disabled = !requests.length;
  }

  // Every notice goes through here, so a storage failure stays on screen.
  function notice(html) {
    $('requestNotice').innerHTML = (storageOk ? '' : storageWarning()) + html;
  }

  function addRequest() {
    var text = tidy($('requestText').value);
    if (!text) { notice(banner('bad', '要望を入力してください。')); return; }
    requests.push({ at: stamp(new Date()), text: text });
    $('requestText').value = '';
    saveRequests();
    notice('');
    renderRequests();
  }

  function legacyCopy(text) {
    var area = document.createElement('textarea');
    area.value = text;
    area.setAttribute('readonly', '');
    area.style.position = 'fixed';
    area.style.opacity = '0';
    document.body.appendChild(area);
    area.select();
    var copied = false;
    try { copied = document.execCommand('copy'); }
    finally { document.body.removeChild(area); }
    return copied;
  }

  // Last resort: put the text where the user can select it by hand.
  function showFallback(boxId, text) {
    $(boxId).innerHTML = banner('warn', '自動でコピーできませんでした。',
      '下の欄を全選択してコピーしてください。')
      + '<textarea readonly rows="8"></textarea>';
    var area = $(boxId).querySelector('textarea');
    area.value = text;
    area.focus();
    area.select();
  }

  /* Clipboard API first, then execCommand, then the text shown for a manual
   * copy in boxId. done() runs only when the text really was copied;
   * failed() runs when it fell through to the manual box. */
  function copyText(text, boxId, done, failed) {
    $(boxId).innerHTML = '';
    function fallback() {
      var copied = false;
      try { copied = legacyCopy(text); }
      catch (error) {
        if (window.console) { console.error('Gibson: copy failed', error); }
      }
      if (copied) { done(); } else { failed(); showFallback(boxId, text); }
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(done, function (error) {
        if (window.console) { console.warn('Gibson: clipboard API refused', error); }
        fallback();
      });
    } else {
      fallback();
    }
  }

  function copyRequests() {
    if (!requests.length) { return; }
    var count = requests.length;
    var pending = tidy($('requestText').value)
      ? banner('warn', '追加されていない文章が入力欄に残っています。') : '';
    copyText(formatRequests(requests, BUILT_AT), 'copyFallback',
      function () { notice(banner('ok', '全 ' + count + ' 件をコピーしました。') + pending); },
      function () { notice(pending); });
  }

  function clearRequests() {
    if (!requests.length) { return; }
    if (!window.confirm(requests.length + ' 件のリクエストをすべて削除します。')) {
      return;
    }
    requests = [];
    saveRequests();
    $('copyFallback').innerHTML = '';
    notice(banner('ok', 'リクエストをすべて削除しました。'));
    renderRequests();
  }

  // ===================================================================
  // Result output
  // ===================================================================

  // The last successful result of each block, with the inputs it was
  // calculated from. A failed calculation clears its block, so an old result
  // never goes into the notebook next to inputs it did not come from.
  var report = { reaction: null, tm: null };

  function hasReport() { return !!(report.reaction || report.tm); }

  function setReport(block, data) {
    report[block] = data;
    var any = hasReport();
    $('copyReport').disabled = !any;
    $('printReport').disabled = !any;
    $('reportNotice').innerHTML = any ? ''
      : '<p class="hint">まだ結果がありません。</p>';
    // A preview of an earlier copy no longer matches the results.
    $('reportPreview').innerHTML = '';
  }

  function showPreview(text) {
    $('reportPreview').innerHTML = '<details open><summary>コピーした内容</summary>'
      + '<pre class="preview"></pre></details>';
    $('reportPreview').querySelector('pre').textContent = text;
  }

  function copyReport() {
    if (!hasReport()) { return; }
    var text = formatReport(report, BUILT_AT, stamp(new Date()));
    $('reportPreview').innerHTML = '';
    copyText(text, 'reportFallback',
      function () {
        $('reportNotice').innerHTML = banner('ok', 'Markdown をコピーしました。');
        showPreview(text);
      },
      // The manual-copy box already shows the whole text.
      function () { $('reportNotice').innerHTML = ''; });
  }

  /* The printed page is the report alone. It is filled on beforeprint, so
   * Cmd+P prints it too; with no result yet the page prints as it is. */
  function preparePrint() {
    var on = hasReport();
    document.body.classList.toggle('print-report', on);
    $('printArea').innerHTML = on
      ? reportToHtml(tablesOnly(formatReport(report, BUILT_AT, stamp(new Date())))) : '';
  }

  // ===================================================================
  // Saved conditions (requested 2026-10-10)
  // ===================================================================

  /* Everything typed into A-2 and B, kept in this browser so a set of
   * conditions can be brought back later. Results are not saved: they are
   * recalculated from the inputs, so they can never go stale. */
  var CONDITION_KEY = 'gibson.conditions.v1';
  var conditions = [];
  var conditionStorageOk = true;

  function conditionWarning() {
    return conditionStorageOk ? '' : banner('warn', 'この画面を閉じると保存した条件は消えます'
      + '（ブラウザに保存できませんでした）。');
  }

  function loadConditions() {
    try {
      var raw = window.localStorage.getItem(CONDITION_KEY);
      if (raw === null) { return; }
      var list = JSON.parse(raw);
      if (!Array.isArray(list)) { throw new Error('saved conditions are not a list'); }
      conditions = list.filter(function (entry) {
        return entry && typeof entry.savedAt === 'string' && Array.isArray(entry.rows);
      });
    } catch (error) {
      conditionStorageOk = false;
      if (window.console) { console.error('Gibson: conditions not loaded', error); }
    }
  }

  function storeConditions() {
    if (!conditionStorageOk) { return; }
    try {
      window.localStorage.setItem(CONDITION_KEY, JSON.stringify(conditions));
    } catch (error) {
      conditionStorageOk = false;
      if (window.console) { console.error('Gibson: conditions not saved', error); }
    }
  }

  function conditionNotice(html) { $('conditionNotice').innerHTML = conditionWarning() + html; }

  function renderConditions() {
    $('conditionList').innerHTML = conditions.length
      ? conditions.map(function (entry, i) {
          return '<div class="request"><div class="when">No.' + (i + 1) + '　'
            + esc(entry.savedAt) + (entry.name ? ' — <strong>' + esc(entry.name) + '</strong>' : '')
            + '</div><div class="body">' + esc(conditionSummary(entry)) + '</div>'
            + '<div class="actions"><button type="button" class="tiny" data-load="' + i
            + '">読み込む</button><button type="button" class="tiny" data-drop="' + i
            + '">削除</button></div></div>';
        }).join('')
      : '<p class="hint">保存した条件はまだありません。</p>';
    $('conditionCount').textContent = '（' + conditions.length + ' 件）';
    Array.prototype.forEach.call($('conditionList').querySelectorAll('[data-load]'),
      function (button) {
        button.addEventListener('click', function () {
          restoreCondition(conditions[parseInt(button.getAttribute('data-load'), 10)]);
        });
      });
    Array.prototype.forEach.call($('conditionList').querySelectorAll('[data-drop]'),
      function (button) {
        button.addEventListener('click', function () {
          var index = parseInt(button.getAttribute('data-drop'), 10);
          var entry = conditions[index];
          if (!window.confirm('No.' + (index + 1) + '（' + entry.savedAt
                              + ' に保存）の条件を削除します。')) { return; }
          conditions.splice(index, 1);
          storeConditions();
          conditionNotice(banner('ok', '削除しました。'));
          renderConditions();
        });
      });
  }

  function saveCondition() {
    conditions.push({ savedAt: stamp(new Date()), name: tidy($('conditionName').value),
                      mode: ratioMode() ? 'ratio' : 'manual',
                      rows: readReaction(),
                      tmText: $('tmText').value });
    $('conditionName').value = '';
    storeConditions();
    // The number is the entry's place in the saved list below.
    conditionNotice(banner('ok', '保存した条件の一覧に No.' + conditions.length
      + ' として追加しました。'));
    renderConditions();
  }

  function restoreCondition(entry) {
    document.querySelector('input[name="amountMode"][value="'
      + (entry.mode === 'ratio' ? 'ratio' : 'manual') + '"]').checked = true;
    $('insertCount').value = String(Math.min(Math.max(entry.rows.length - 1, 1), MAX_INSERTS));
    buildReactionRows();
    Array.prototype.forEach.call($('reactionRows').rows, function (tr, i) {
      REACTION_FIELDS.forEach(function (field) {
        var input = tr.querySelector('[data-field="' + field + '"]');
        if (input) { input.value = (entry.rows[i] && entry.rows[i][field]) || ''; }
      });
    });
    $('tmText').value = entry.tmText || '';
    $('tmResult').innerHTML = '';
    setReport('tm', null);
    applyMode();
    conditionNotice(banner('ok', 'No.' + (conditions.indexOf(entry) + 1) + '（'
      + entry.savedAt + ' に保存）の条件を読み込みました。「計算・判定」を押すと結果が出ます。'));
  }

  // ===================================================================
  // Self test
  // ===================================================================

  function runSelfTests() {
    var passed = 0, failures = {};
    PBTESTS.cases.forEach(function (item) {
      try { item.fn(); passed++; }
      catch (error) { failures[item.name] = error.message; }
    });
    var failed = PBTESTS.cases.length - passed;
    var html = failed
      ? banner('bad', passed + ' / ' + PBTESTS.cases.length + ' 通過（'
               + failed + ' 件失敗）', 'このツールの計算は信用できません。')
      : banner('ok', passed + ' / ' + PBTESTS.cases.length + ' すべて通過',
               '量の換算・Tm・リクエストの書式が、手計算した期待値と一致しました。');
    html += '<div class="scroll"><table><tr><th></th><th>項目</th><th>結果</th></tr>';
    PBTESTS.cases.forEach(function (item) {
      var failure = failures[item.name];
      html += '<tr><td>' + (failure ? '❌' : '✅') + '</td><td>' + esc(item.name)
        + '</td><td>' + (failure ? esc(failure) : 'OK') + '</td></tr>';
    });
    $('testResult').innerHTML = html + '</table></div>';
  }

  // ===================================================================
  // Wiring
  // ===================================================================

  $('mixLabel').textContent = 'G.A. M ×' + GIBSON.MASTER_MIX_FOLD;
  $('totalVolume').textContent = GIBSON.TOTAL_UL;
  for (var n = 1; n <= MAX_INSERTS; n++) {
    $('insertCount').insertAdjacentHTML('beforeend',
      '<option value="' + n + '">' + n + ' 本</option>');
  }
  // ---- Tab bar follows the scroll (requested 2026-10-12) ---------------
  var navLinks = Array.prototype.slice.call(document.querySelectorAll('#sectionNav a'));
  var navTargets = navLinks.map(function (a) {
    return document.getElementById(a.getAttribute('href').slice(1));
  });
  var navPending = false;
  // A tab that was clicked stays highlighted until the user scrolls by
  // hand: the short sections near the end can never scroll up to the bar,
  // so without this the click would land on "自己テスト" instead.
  var navLocked = -1;
  function highlightSection() {
    navPending = false;
    var bar = document.querySelector('.appbar').offsetHeight;
    // Sections stop 16 px below the bar when jumped to (scroll-margin-top),
    // and the line sits just under that, so a jumped-to section counts.
    document.documentElement.style.setProperty('--bar-h', bar + 'px');
    var doc = document.documentElement;
    var atBottom = window.innerHeight + window.scrollY >= doc.scrollHeight - 2;
    var current = navLocked >= 0 ? navLocked
      : sectionInView(navTargets.map(function (el) {
          return el.getBoundingClientRect().top;
        }), bar + 24, atBottom);
    navLinks.forEach(function (a, i) { a.classList.toggle('active', i === current); });
  }
  navLinks.forEach(function (a, i) {
    a.addEventListener('click', function () { navLocked = i; highlightSection(); });
  });
  ['wheel', 'touchmove', 'keydown'].forEach(function (type) {
    window.addEventListener(type, function () { navLocked = -1; }, { passive: true });
  });
  // One update per frame, however fast the scroll events come.
  window.addEventListener('scroll', function () {
    if (!navPending) { navPending = true; window.requestAnimationFrame(highlightSection); }
  }, { passive: true });
  window.addEventListener('resize', highlightSection);
  highlightSection();

  $('guideBox').innerHTML = guideHtml();
  buildReactionRows();
  $('insertCount').addEventListener('change', buildReactionRows);
  Array.prototype.forEach.call(document.querySelectorAll('input[name="amountMode"]'),
    function (radio) { radio.addEventListener('change', applyMode); });
  // Once a number box is left, show what will be read: "．８" becomes "0.8".
  function tidyNumberBox(event) {
    var box = event.target;
    // Names are text, not numbers: leave "１番" alone.
    if (box.tagName === 'INPUT' && box.type === 'text' && !box.readOnly
        && box.getAttribute('data-field') !== 'name') {
      box.value = GIBSON.normalizeNumberText(box.value);
    }
  }
  $('reactionTable').addEventListener('change', tidyNumberBox);
  // Editing any input invalidates the numbers next to it.
  $('reactionTable').addEventListener('input', clearReaction);
  $('pCalc').addEventListener('click', calcReaction);
  $('recButton').addEventListener('click', recommendAmounts);
  $('tmCalc').addEventListener('click', calcTm);
  $('addRequest').addEventListener('click', addRequest);
  $('copyRequests').addEventListener('click', copyRequests);
  $('clearRequests').addEventListener('click', clearRequests);
  $('testButton').addEventListener('click', runSelfTests);
  $('saveCondition').addEventListener('click', saveCondition);
  $('copyReport').addEventListener('click', copyReport);
  // Filled here as well as on beforeprint: not every browser fires that
  // event before laying out the page.
  $('printReport').addEventListener('click', function () {
    preparePrint();
    window.print();
  });
  window.addEventListener('beforeprint', preparePrint);
  window.addEventListener('afterprint', function () {
    document.body.classList.remove('print-report');
  });
  setReport('reaction', null);
  loadConditions();
  conditionNotice('');
  renderConditions();
  loadRequests();
  notice('');
  renderRequests();

  return exported;
})();
