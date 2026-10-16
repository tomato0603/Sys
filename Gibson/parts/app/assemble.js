/* ===================================================================
 * Building the new region out of named pieces.
 *
 * This unit is what separates Construct Builder from Primer Builder. Both
 * are built from the same parts; only this one is added. Everything else --
 * the engine, the sequence view, the reading frame, the primer search, the
 * check, the files -- is the same code, so the two pages cannot drift.
 *
 * What it adds is the step in the middle: instead of typing the peptide to
 * put in, assemble it from the pieces already imported from Protein Builder
 * (tags, linkers, domains), in order, and hand the result to the field the
 * design reads. The pieces are the same stored list the colour marks use.
 *
 * It owns nothing else. The moment the sequence is handed over it is an
 * ordinary amino-acid replacement, checked for reading frame and verified
 * like any other.
 * =================================================================== */
var PBASM = (function () {
  'use strict';

  var AMINO_ONLY = /^[ACDEFGHIKLMNPQRSTVWY]+$/;

  // Pieces in the order they will be joined. Not saved anywhere: a build in
  // progress is a thought, and the result goes into the design.
  var pieces = [];

  function $(id) { return document.getElementById(id); }
  function esc(t) {
    return String(t === undefined || t === null ? '' : t)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }
  function clean(text) {
    return String(text || '').toUpperCase().replace(/[^A-Z]/g, '');
  }
  function joined() {
    return pieces.map(function (p) { return p.sequence; }).join('');
  }

  function move(index, by) {
    var to = index + by;
    if (to < 0 || to >= pieces.length) { return; }
    var moved = pieces.splice(index, 1)[0];
    pieces.splice(to, 0, moved);
    render();
  }

  function addFromLibrary() {
    var chooser = $('asmPick');
    var entry = (window.PBUI.proteins() || [])[parseInt(chooser.value, 10)];
    if (!entry) { return; }
    pieces.push({ name: entry.name, sequence: entry.sequence });
    render();
  }

  function addTyped() {
    var field = $('asmTyped'), sequence = clean(field.value);
    if (!sequence) { return; }
    if (!AMINO_ONLY.test(sequence)) {
      $('asmNote').textContent = 'アミノ酸として読めない文字があります。'
        + '1 文字表記で入力してください。';
      return;
    }
    pieces.push({ name: ($('asmName').value || '').trim() || '手入力',
                  sequence: sequence });
    field.value = '';
    $('asmName').value = '';
    render();
  }

  function rows() {
    if (!pieces.length) {
      return '<p class="hint">まだ何も並んでいません。'
        + '下から部品を足してください。</p>';
    }
    var at = 0;
    var html = '<div class="scroll"><table><tr><th>#</th><th>部品</th>'
      + '<th>残基</th><th>位置</th><th></th></tr>';
    pieces.forEach(function (piece, index) {
      var from = at + 1;
      at += piece.sequence.length;
      html += '<tr><td>' + (index + 1) + '</td><td>' + esc(piece.name)
        + '</td><td>' + piece.sequence.length + '</td><td>' + from + '..' + at
        + '</td><td>'
        + '<button type="button" class="tiny" data-asmup="' + index
        + '"' + (index === 0 ? ' disabled' : '') + '>↑</button> '
        + '<button type="button" class="tiny" data-asmdown="' + index + '"'
        + (index === pieces.length - 1 ? ' disabled' : '') + '>↓</button> '
        + '<button type="button" class="tiny" data-asmdel="' + index
        + '">外す</button></td></tr>';
    });
    return html + '</table></div>';
  }

  function render() {
    var box = $('assembleStep');
    if (!box) { return; }
    var library = window.PBUI.proteins() || [];
    var sequence = joined();

    var html = '<h2>部品を並べて配列を作る（任意）</h2>'
      + '<p class="hint">Protein Builder から取り込んだ部品を並べます。'
      + 'できた配列は下の「新しい配列」に入れて、そのまま設計できます。</p>'
      + rows();

    if (sequence) {
      html += '<table><tr><th>できた配列</th><td class="seq">'
        + sequence.length + ' 残基 ' + esc(sequence) + '</td></tr></table>'
        + '<div class="actions">'
        + '<button type="button" class="btn" id="asmUse">'
        + 'この配列を「新しい配列」に入れる</button>'
        + '<button type="button" class="tiny" id="asmCopy">コピー</button>'
        + '<button type="button" class="tiny" id="asmClear">全部外す</button>'
        + '</div>';
    }

    html += '<h3>部品を足す</h3>';
    if (library.length) {
      html += '<div class="row"><div>'
        + '<label for="asmPick">取り込んだ配列から</label><select id="asmPick">';
      library.forEach(function (entry, index) {
        html += '<option value="' + index + '">' + esc(entry.name) + '（'
          + entry.sequence.length + ' 残基）</option>';
      });
      html += '</select></div><div><label>&nbsp;</label>'
        + '<button type="button" id="asmAdd">末尾に足す</button></div></div>';
    } else {
      html += '<p class="hint">取り込んだ配列がまだありません。'
        + '上の「Protein Builder のタンパク質配列」で '
        + '<code>.json</code> を読み込んでください。</p>';
    }
    html += '<div class="row"><div>'
      + '<label for="asmName">名前<small>省略できます</small></label>'
      + '<input type="text" id="asmName" placeholder="GS linker"></div>'
      + '<div><label for="asmTyped">直接入力'
      + '<small>アミノ酸 1 文字表記</small></label>'
      + '<input type="text" id="asmTyped" placeholder="GGGGS"></div>'
      + '<div><label>&nbsp;</label>'
      + '<button type="button" id="asmAddTyped">末尾に足す</button></div></div>'
      + '<p class="hint" id="asmNote"></p>';

    box.innerHTML = html;
    wire();
  }

  function wire() {
    var box = $('assembleStep');
    [['asmAdd', addFromLibrary], ['asmAddTyped', addTyped],
     ['asmClear', function () { pieces = []; render(); }],
     ['asmCopy', function () { window.PBUI.copy(joined(), $('asmCopy')); }],
     ['asmUse', function () { window.PBUI.usePeptide(joined()); }]
    ].forEach(function (pair) {
      var button = $(pair[0]);
      if (button) { button.addEventListener('click', pair[1]); }
    });
    [['data-asmup', -1], ['data-asmdown', 1]].forEach(function (pair) {
      Array.prototype.forEach.call(box.querySelectorAll('[' + pair[0] + ']'),
        function (button) {
          button.addEventListener('click', function () {
            move(parseInt(button.getAttribute(pair[0]), 10), pair[1]);
          });
        });
    });
    Array.prototype.forEach.call(box.querySelectorAll('[data-asmdel]'),
      function (button) {
        button.addEventListener('click', function () {
          pieces.splice(parseInt(button.getAttribute('data-asmdel'), 10), 1);
          render();
        });
      });
  }

  document.addEventListener('DOMContentLoaded', function () {
    // The two builds must be tellable apart at a glance, and in the browser's
    // tab and history.
    document.title = 'Construct Builder';
    var heading = document.querySelector('h1');
    if (heading) { heading.textContent = 'Construct Builder'; }
    render();
  });

  return { refresh: render, pieces: function () { return pieces.slice(); } };
})();
if (typeof module !== 'undefined' && module.exports) { module.exports = PBASM; }
