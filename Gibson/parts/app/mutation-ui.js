/* Primer Builder: the mutation-notation mode.
 *
 * The change is stated as K42A / ins@50:GGSGGS / del45-50 / rep45-50:GS, or as
 * a whole target protein to diff against. Several edits can sit too far apart
 * to reach with one primer pair, so this is the mode that splits the work into
 * rounds. The other mode -- picking a range on screen -- lives in ui.js.
 *
 * The plasmid arrives through setRecord(); this unit never reads a file
 * itself. Copying and saving are the shared ones, reached via window.PBUI.
 */
var PBMUT = (function () {
  'use strict';

  var state = { record: null, feature: null, protein: '', run: null };

  function esc(text) {
    return String(text === undefined || text === null ? '' : text)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }
  function $(id) { return document.getElementById(id); }
  function show(id, visible) {
    $(id).style.display = visible ? '' : 'none';
  }
  function banner(kind, title, body) {
    return '<div class="banner ' + kind + '"><strong>' + esc(title)
      + '</strong>' + (body || '') + '</div>';
  }
  function message(text, kind) {
    $('messages').innerHTML = text
      ? banner(kind || 'bad', kind === 'ok' ? 'できました' : '入力を確認してください',
               esc(text))
      : '';
  }
  function copyButton(text, label) {
    return '<button type="button" class="tiny copy" data-copy="' + esc(text)
      + '">' + esc(label || 'コピー') + '</button>';
  }

  // ---- step 2: choose the CDS ----------------------------------------

  function renderCdsChooser() {
    var features = state.record.cdsFeatures();
    var html = '';
    features.forEach(function (feature, index) {
      html += '<label><input type="radio" name="cds" value="' + index + '"'
        + (index === 0 ? ' checked' : '') + '> ' + esc(feature.label())
        + ' <small>' + esc(feature.locationText()) + '</small></label>';
    });
    html += '<label><input type="radio" name="cds" value="manual"'
      + (features.length ? '' : ' checked') + '> 範囲を自分で指定する</label>'
      + '<input type="text" id="cdsManual" placeholder="例: 101..700 / '
      + 'complement(100..200) / join(5700..5715,1..100)">';
    $('cdsChooser').innerHTML = html;
    Array.prototype.forEach.call(
      document.querySelectorAll('input[name=cds]'), function (input) {
        input.addEventListener('change', applyCds);
      });
    $('cdsManual').addEventListener('change', function () {
      var manual = document.querySelector('input[name=cds][value=manual]');
      if (manual) { manual.checked = true; }
      applyCds();
    });
    applyCds();
  }

  function applyCds() {
    var checked = document.querySelector('input[name=cds]:checked');
    var choice = checked ? checked.value : '';
    var manual = ($('cdsManual').value || '').trim();
    try {
      if (choice === 'manual') {
        if (!manual) {
          state.feature = null;
          $('proteinView').innerHTML = banner('warn',
            'CDS の範囲を入力してください',
            (state.record.cdsFeatures().length
              ? '「範囲を自分で指定する」を選んだので、上の欄に範囲を入力してください。'
              : 'このファイルには CDS 注釈がありません。改変対象の遺伝子の範囲を'
                + '上の欄に入力してください。')
            + '<br>書き方: <code>101..700</code>（順方向）／'
            + '<code>complement(100..200)</code>（逆方向）／'
            + '<code>join(5700..5715,1..100)</code>（原点をまたぐ）'
            + '<br>範囲が分からない場合は、ApE でその遺伝子の開始位置と終了位置を'
            + '確認してください。');
          $('cdsManual').focus();
          return;
        }
        var parsed = PB.parseLocation(manual);
        state.feature = new PB.Feature('CDS', parsed.segments, parsed.strand,
                                       { label: '画面で指定' });
      } else {
        state.feature = state.record.cdsFeatures()[parseInt(choice, 10)];
      }
      var translated = PB.translate(state.record.featureSequence(state.feature));
      state.protein = translated.protein;
      var warning = translated.trailing
        ? '<p class="hint">注意: CDS の長さが 3 の倍数ではありません（余り '
          + translated.trailing + ' 塩基）。範囲を確認してください。</p>' : '';
      $('proteinView').innerHTML = '<h3>この CDS の翻訳（残基番号つき）</h3>'
        + '<p class="hint">変異表記の残基番号はこの番号です。'
        + '残基 1 は開始メチオニンです。</p>' + warning
        + '<div class="protein">' + numberedProtein(state.protein) + '</div>';
      message('');
      show('mnEdit', true);
    } catch (error) {
      state.feature = null;
      $('proteinView').innerHTML = '';
      message(error.message);
    }
  }

  function numberedProtein(protein) {
    var body = protein.replace(/\*+$/, '');
    var rows = [];
    for (var start = 0; start < body.length; start += 30) {
      var chunk = body.substr(start, 30), spaced = '';
      for (var i = 0; i < chunk.length; i += 10) {
        spaced += (i ? ' ' : '') + chunk.substr(i, 10);
      }
      var label = String(start + 1);
      while (label.length < 4) { label = ' ' + label; }
      rows.push('<span class="pos">' + label + '</span>  ' + esc(spaced));
    }
    return rows.join('\n');
  }

  // ---- step 3: design -------------------------------------------------

  function integerField(id, fallback) {
    var text = ($(id).value || '').trim();
    if (!text) { return fallback; }
    if (!/^\d+$/.test(text)) {
      throw new Error($(id).getAttribute('data-label')
                      + ' は整数で入力してください: ' + text);
    }
    var value = parseInt(text, 10);
    if (value <= 0) {
      throw new Error($(id).getAttribute('data-label') + ' は 1 以上で入力してください');
    }
    return value;
  }

  /* The usage floor, as a fraction. Read separately from integerField
   * because 0 is a legitimate value here -- it means "no band at all" -- and
   * integerField rejects anything below 1. */
  function usageFloor() {
    var raw = ($('minUsage').value || '').trim();
    if (!raw) { return 0.8; }
    if (!/^\d+$/.test(raw)) {
      throw new Error('使用頻度の下限は 0 から 100 の整数で入力してください: ' + raw);
    }
    var value = parseInt(raw, 10);
    if (value > 100) {
      throw new Error('使用頻度の下限は 100 以下で入力してください');
    }
    return value / 100;
  }

  function design(diagnosticEnzyme) {
    if (!state.feature) {
      message('改変対象の CDS が決まっていません。'
              + '上の「3. 改変対象の遺伝子（CDS）」で範囲を指定してください。');
      $('mnCds').scrollIntoView({ behavior: 'smooth', block: 'start' });
      return;
    }
    var mutations = ($('mutations').value || '').trim();
    var target = ($('target').value || '').trim();
    if (target) {
      var body = target.charAt(0) === '>'
        ? target.slice(target.indexOf('\n') + 1) : target;
      target = body.replace(/[\s\d]/g, '').toUpperCase();
    }
    if (!!mutations === !!target) {
      message('変異表記か目的タンパク質の配列を、どちらか一方だけ入力してください。');
      return;
    }
    try {
      state.run = PB.designRun(state.record, state.feature, {
        mutations: mutations || null, targetProtein: target || null,
        maxPrimerLen: integerField('maxPrimerLen', PB.DEFAULT_MAX_PRIMER_LEN),
        minAnnealLen: integerField('minAnnealLen', PB.DEFAULT_MIN_ANNEAL_LEN),
        maxAnnealLen: integerField('maxAnnealLen', PB.DEFAULT_MAX_ANNEAL_LEN),
        minUsage: usageFloor(),
        diagnosticEnzyme: diagnosticEnzyme
          || (($('diagnostic').value || '').trim() || null)
      });
    } catch (error) {
      state.run = null;
      $('result').innerHTML = '';
      message(error.message);
      return;
    }
    message('');
    if (diagnosticEnzyme) { $('diagnostic').value = diagnosticEnzyme; }
    renderResult();
    $('result').scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  function renderResult() {
    var run = state.run, strip = function (p) { return p.replace(/\*+$/, ''); };
    var html = '';
    html += run.ok
      ? banner('ok', 'すべてのラウンドで設計と検算が通りました',
               '下のプライマーを発注してください。')
      : banner('bad', '問題の残るラウンドがあります',
               '各ラウンドの内容を確認してください。');
    run.notes.forEach(function (note) {
      html += banner('warn', '注意', esc(note));
    });

    html += '<h2>編集内容</h2><table><tr><th>#</th><th>指定</th><th>内容</th></tr>';
    run.edits.forEach(function (edit, index) {
      html += '<tr><td>' + (index + 1) + '</td><td class="mono">'
        + esc(edit.source || '差分') + '</td><td>' + esc(edit.describe())
        + '</td></tr>';
    });
    html += '</table><p>タンパク質 ' + strip(run.protein).length + ' 残基 → '
      + strip(run.newProtein).length + ' 残基。逆PCR ' + run.rounds.length
      + ' ラウンド（5′ テールの上限 ' + run.budget + ' nt）。</p>';

    html += '<h2>プライマー</h2>'
      + '<div class="scroll"><table><tr><th>名前</th><th>配列 (5′→3′)</th>'
      + '<th>長さ</th><th>Tm</th><th>GC</th><th></th></tr>';
    PB.primerRows(run).forEach(function (row) {
      html += '<tr><td class="mono">' + esc(row.name) + '</td><td class="seq">'
        + esc(row.sequence) + '</td><td>' + row.length + ' nt</td><td>'
        + row.tm + '℃</td><td>' + row.gc.toFixed(0) + '%</td><td>必要</td><td>'
        + copyButton(row.sequence) + '</td></tr>';
    });
    html += '</table></div>';
    html += '<p>' + copyButton(PB.primerRows(run).map(function (r) {
      return r.name + '\t' + r.sequence;
    }).join('\n'), '全部まとめてコピー（名前と配列、タブ区切り）') + '</p>';

    html += '<h2>ダウンロード</h2><p>'
      + '<button type="button" class="btn" data-download="ape">'
      + '改変後のプラスミド (.ape)</button> '
      + '<button type="button" class="btn" data-download="md">'
      + '設計レポート (.md)</button>'
      + (run.ok ? ' <button type="button" class="btn" id="saveMutationDesign">'
                  + 'この設計を保存する</button>' : '')
      + '</p><p class="hint">.ape は ApE でそのまま開けます。'
      + '親プラスミドの COMMENT を保ったうえで、今回の改変内容とプライマー配列が'
      + ' COMMENT として追記されます。<br>'
      + 'ブラウザがダウンロードを止める場合は、下の「テキストで開く」から'
      + 'コピーしてください。</p><p>'
      + '<button type="button" class="tiny" data-showtext="ape">'
      + '.ape をテキストで開く</button> '
      + '<button type="button" class="tiny" data-showtext="md">'
      + 'レポートをテキストで開く</button></p>';
      // The panel itself lives in head.html, below #result. There used to be
      // a second one with the same id here, so $('textPanel') picked
      // whichever came first in the document (2026-09-27).

    html += '<h2>ラウンドごとの詳細</h2>';
    run.rounds.forEach(function (round) { html += roundCard(run, round); });
    $('result').innerHTML = html;
    wireResultButtons();
  }

  function roundCard(run, round) {
    var html = '<div class="card"><h3>ラウンド ' + round.index + ' &mdash; '
      + esc(round.group.edits.map(function (e) {
          return e.source || '差分';
        }).join(', ')) + '</h3>';

    if (!round.pair) {
      round.problems.forEach(function (problem) {
        html += banner('bad', '設計できません', esc(problem));
      });
      return html + '</div>';
    }

    var pair = round.pair, site = pair.site, v = round.verification;
    var strip = function (p) { return p.replace(/\*+$/, ''); };
    if (v && v.ok) {
      html += banner('ok', '検算に通りました',
        'プライマーの配列から環状化後のプラスミドを組み立て直し、'
        + 'CDS を配列検索で特定して翻訳したところ、目的のタンパク質（'
        + strip(v.observedProtein).length + ' 残基）と完全に一致しました。');
    } else {
      var problems = round.problems.concat(v ? v.problems : []);
      html += banner('bad', '問題があります',
                     problems.map(esc).join('<br>'));
    }
    if (v) {
      v.notes.forEach(function (note) {
        html += banner('warn', '注意', esc(note));
      });
    }

    html += '<table>'
      + '<tr><th>除去する領域</th><td>プラスミド ' + site.leftPos
      + ' の直後から ' + site.rightPos + ' の直前（' + site.removedLength
      + ' bp）</td></tr>'
      + '<tr><th>挿入する配列</th><td class="seq">' + site.insertLength + ' bp '
      + esc(site.insert) + ' ' + copyButton(site.insert) + '</td></tr>'
      + '<tr><th>テールの分配</th><td>R に ' + pair.split + ' nt / F に '
      + (site.insertLength - pair.split) + ' nt</td></tr>'
      + '<tr><th>予想産物長</th><td>' + pair.productLength + ' bp</td></tr>'
      + '</table>';

    html += '<h3>プライマー</h3><div class="scroll"><table>'
      + '<tr><th></th><th>配列 (5′→3′)</th><th>長さ</th><th>5′ テール</th>'
      + "<th>3′ アニール域</th><th>Tm</th><th>GC</th><th>3′ 末端</th>"
      + '<th>結合部位</th><th></th></tr>';
    [[round.forwardName, pair.forward], [round.reverseName, pair.reverse]]
      .forEach(function (entry) {
        var primer = entry[1];
        html += '<tr><td class="mono">' + esc(entry[0]) + '</td>'
          + '<td class="seq">' + esc(primer.sequence) + '</td><td>'
          + primer.length + ' nt</td><td>' + primer.tail.length + ' nt</td><td>'
          + primer.anneal.length + ' nt</td><td>' + primer.tm + '℃</td><td>'
          + primer.gcPercent.toFixed(0) + '%</td><td>' + esc(primer.terminalBase)
          + '</td><td>' + primer.bindingSites + ' 箇所'
          + (primer.unique ? '' : ' ⚠') + '</td><td>'
          + copyButton(primer.sequence) + '</td></tr>';
      });
    html += '</table></div><p class="hint">Tm 差 ' + pair.tmDifference
      + '℃。3′ アニール域は鋳型と完全一致で、編集はすべて 5′ テールに載っています。'
      + '「結合部位 1 箇所」はプラスミド全体（両鎖・原点跨ぎ込み）で一意である'
      + 'ことの確認です。</p>';

    html += PBUI.pcrForm();

    if (round.diagnostics.length) {
      html += '<h3>診断用制限酵素サイト</h3><p class="hint">'
        + '親プラスミドに無く、改変後に 1 箇所だけ現れる酵素です。'
        + 'アミノ酸配列は変わりません。採用すると挿入配列が書き換わり、'
        + 'プライマーも再設計されます。</p><table><tr><th>酵素</th>'
        + '<th>認識配列</th><th>変える塩基</th><th></th></tr>';
      round.diagnostics.forEach(function (candidate) {
        var adopted = round.adoptedDiagnostic
          && round.adoptedDiagnostic.enzyme === candidate.enzyme;
        html += '<tr><td>' + esc(candidate.enzyme) + '</td><td class="mono">'
          + esc(candidate.recognition) + '</td><td>' + candidate.changedBases
          + ' 塩基</td><td>' + (adopted ? '<strong>採用中</strong>'
            : '<button type="button" class="tiny" data-adopt="'
              + esc(candidate.enzyme) + '">これを採用して再設計</button>')
          + '</td></tr>';
      });
      html += '</table>';
    }

    if (round.alternatives.length) {
      html += '<details><summary>代替候補 ' + round.alternatives.length
        + ' 件（スコアは減点方式、小さいほど良い）</summary>'
        + '<div class="scroll"><table><tr><th></th><th>F 長</th><th>F Tm</th>'
        + '<th>R 長</th><th>R Tm</th><th>Tm 差</th><th>スコア</th></tr>';
      [['採用', pair]].concat(round.alternatives.map(function (a) {
        return ['', a];
      })).forEach(function (entry) {
        var candidate = entry[1];
        html += '<tr><td>' + entry[0] + '</td><td>' + candidate.forward.length
          + '</td><td>' + candidate.forward.tm + '℃</td><td>'
          + candidate.reverse.length + '</td><td>' + candidate.reverse.tm
          + '℃</td><td>' + candidate.tmDifference + '℃</td><td>'
          + candidate.score.toFixed(1) + '</td></tr>';
      });
      html += '</table></div><p class="hint">候補は 3′ アニーリング領域の'
        + '長さの組み合わせで、テールの内容は同じです。</p></details>';
    }
    return html + '</div>';
  }

  // ---- downloads ------------------------------------------------------

  /* The date goes in front of the FILE name only (asked for 2026-09-25):
   * one download folder holds many designs. The record's own name stays as
   * the construct name, because that is what ApE shows on the LOCUS line and
   * what the next design will be built from. */
  function fileBaseName() {
    var run = state.run;
    return PBUI.stamp() + '_'
      + ((run && run.constructName) || 'design').replace(/[^\w.\-]/g, '_');
  }

  function outputFiles() {
    var run = state.run;
    if (!run) { return {}; }
    // The date in front: one download folder, many designs (2026-09-25).
    var base = (run.constructName || 'design').replace(/[^\w.\-]/g, '_');
    var fileBase = fileBaseName();
    var files = { md: { name: fileBase + '_design.md',
                        text: PB.designMarkdown(run) } };
    for (var i = 0; i < run.rounds.length; i++) {
      var round = run.rounds[i];
      if (!round.verification) { continue; }
      var record = round.verification.newRecord;
      var suffix = run.rounds.length === 1 ? '' : '_round' + round.index;
      record.name = base + suffix;
      record.headerLines = run.record.headerLines
        .concat(PB.editCommentLines(run, round));
      files[i === 0 ? 'ape' : 'ape' + round.index] = {
        name: fileBase + suffix + '.ape', text: PB.formatGenbank(record) };
      if (run.rounds.length === 1) { break; }
    }
    return files;
  }

  function download(kind) {
    var file = outputFiles()[kind];
    if (!file) { message('ダウンロードできるファイルがありません。'); return; }
    try {
      var blob = new Blob([file.text], { type: 'text/plain;charset=utf-8' });
      var url = URL.createObjectURL(blob);
      var link = document.createElement('a');
      link.href = url; link.download = file.name;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
    } catch (error) {
      showText(kind);
      message('ブラウザがダウンロードを止めました。下のテキストからコピーしてください。');
    }
  }

  function showText(kind) {
    var file = outputFiles()[kind];
    if (!file) { return; }
    $('textPanel').innerHTML = '<h3>' + esc(file.name) + '</h3>'
      + '<p>' + copyButton(file.text, '全文をコピー') + '</p>'
      + '<textarea rows="14" readonly>' + esc(file.text) + '</textarea>';
  }

  function wireResultButtons() {
    var save = document.getElementById('saveMutationDesign');
    if (save) {
      save.addEventListener('click', function () {
        var run = state.run, files = outputFiles();
        PBUI.saveDesign({
          name: run.constructName || run.record.name,
          base: fileBaseName(),
          summary: run.edits.map(function (edit) {
            return edit.source || '差分';
          }).join(', '),
          report: files.md ? files.md.text : '',
          template: PB.formatGenbank(run.record),
          edited: files.ape ? files.ape.text : ''
        });
      });
    }
    Array.prototype.forEach.call(
      $('result').querySelectorAll('[data-copy]'), function (button) {
        button.addEventListener('click', function () {
          PBUI.copy(button.getAttribute('data-copy'), button);
        });
      });
    Array.prototype.forEach.call(
      $('result').querySelectorAll('[data-download]'), function (button) {
        button.addEventListener('click', function () {
          download(button.getAttribute('data-download'));
        });
      });
    Array.prototype.forEach.call(
      $('result').querySelectorAll('[data-showtext]'), function (button) {
        button.addEventListener('click', function () {
          showText(button.getAttribute('data-showtext'));
        });
      });
    Array.prototype.forEach.call(
      $('result').querySelectorAll('[data-adopt]'), function (button) {
        button.addEventListener('click', function () {
          design(button.getAttribute('data-adopt'));
        });
      });
  }

  // ---- entry points used by the shared UI -----------------------------

  /* Called by ui.js whenever a plasmid is loaded. The CDS chooser is rebuilt
   * from scratch: the previous file's features must not survive. */
  function setRecord(record) {
    state.record = record;
    state.feature = null;
    state.protein = '';
    state.run = null;
    renderCdsChooser();
  }

  function clearResult() {
    state.run = null;
    $('result').innerHTML = '';
    $('textPanel').innerHTML = '';
  }

  /* What the request form reports for this mode. How the edit was stated, not
   * the sequence itself. */
  function settings() {
    var parts = [];
    if (state.feature) { parts.push('CDS: ' + state.feature.locationText()); }
    var mutations = ($('mutations').value || '').trim();
    if (mutations) { parts.push('変異表記: ' + mutations); }
    if (($('target').value || '').trim()) { parts.push('目的タンパク質: 入力あり'); }
    if (state.run) { parts.push('ラウンド数: ' + state.run.rounds.length); }
    return parts;
  }

  return { setRecord: setRecord, clearResult: clearResult,
           settings: settings, design: design };
})();
if (typeof module !== 'undefined' && module.exports) { module.exports = PBMUT; }
