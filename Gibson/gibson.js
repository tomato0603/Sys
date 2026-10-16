/* ===================================================================
 * Gibson core -- fragment amounts and primer Tm.
 *
 * Two calculations, and nothing else: converting between ng, pmol and
 * microlitres for the fragments going into a Gibson reaction, and the Tm of a
 * pasted oligo.
 *
 * It sits on top of PB (parts/shared/engine.js) and does NOT carry its own
 * Tm formula. Having the Takara formula in two files is exactly the failure
 * this project has already lived through once (see CLAUDE.md 2.2), and a Tm
 * that disagrees with the one Primer Builder prints for the same oligo would
 * be worse than no Tm at all.
 *
 * The dependency runs one way. PB knows nothing about GIBSON, so the two
 * existing pages are unaffected. PB is not touched while this file loads,
 * only when a function is called, so the load order is free.
 *
 * What is deliberately absent: target pmol, molar ratios, extension times,
 * overlap design. Those are decisions about how to run the experiment, and
 * they belong to the user (CLAUDE.md 0). This file converts; it does not
 * recommend.
 * =================================================================== */
var GIBSON = (function () {
  'use strict';

  // Conventional average molecular weight of one base pair of dsDNA. It is an
  // average, so a GC-skewed fragment is slightly off; the tool claims no more
  // precision than that.
  var DS_DNA_G_PER_MOL_PER_BP = 650;

  // Below this, pipetting error swamps the number we just computed.
  var MIN_PIPETTE_UL = 0.5;

  // Gibson Assembly Master Mix is supplied at 2x, so it is half the reaction.
  var MASTER_MIX_FOLD = 2;
  // The reaction is always 20 uL in the user's lab (given 2026-10-10), so
  // it is not an input.
  var TOTAL_UL = 20;

  /* The lab's own guide for checking hand-entered amounts (given by the user
   * on 2026-10-10). They are checks, not answers: reaction() reports whether
   * the amounts meet them and never picks an amount itself. */
  var VECTOR_NG = [50, 100];
  var INSERT_RATIO_LONG = [2, 3];       // insert ng / vector ng, insert >= 200 bp
  var INSERT_RATIO_SHORT_MIN = 5;       // the same, insert < 200 bp
  var SHORT_INSERT_BP = 200;
  // The smallest volume the user's pipette can take (given 2026-10-10).
  var MIN_AMOUNT_UL = 0.22;
  // Total pmol by number of fragments = vector + inserts, so there is never
  // a 1-fragment reaction. The user corrected the guide's "1-2" to "2-3" on
  // 2026-10-10, which also gives 3 fragments a range.
  var TOTAL_PMOL_RANGES = [
    { min: 2, max: 3, pmol: [0.02, 0.5] },
    { min: 4, max: 6, pmol: [0.2, 1.0] }
  ];

  var PMOL_DIGITS = 3;      // significant figures
  var NG_DIGITS = 4;        // significant figures -- 162.5 ng stays checkable
                            // by hand; 3 digits would print 163.
  var UL_DECIMALS = 2;
  var GC_DECIMALS = 1;

  function GibsonError(message) {
    var error = new Error(message);
    error.name = 'GibsonError';
    return error;
  }

  function engine() {
    if (typeof PB === 'undefined') {
      throw GibsonError('計算エンジン (PB) が読み込まれていません。');
    }
    return PB;
  }

  // =====================================================================
  // Numbers
  // =====================================================================

  /* Number text as typed with a Japanese IME: full-width digits, and "."
   * pressed in Japanese mode, which comes out as "。" or "．". ".8" is also
   * written out as "0.8" (requested 2026-10-10). Anything else is left as
   * it is, for positive() to accept or refuse. */
  function normalizeNumberText(text) {
    var t = String(text === undefined || text === null ? '' : text)
      .replace(/[０-９]/g, function (d) { return String.fromCharCode(d.charCodeAt(0) - 0xFEE0); })
      .replace(/[．。｡]/g, '.').replace(/^\s+|\s+$/g, '');
    return /^\.\d/.test(t) ? '0' + t : t;
  }

  function positive(value, what) {
    if (typeof value === 'string') { value = normalizeNumberText(value); }
    var n = typeof value === 'number' ? value : parseFloat(value);
    if (value === '' || value === null || value === undefined || isNaN(n)) {
      throw GibsonError(what + 'を入力してください。');
    }
    if (!isFinite(n) || n <= 0) {
      throw GibsonError(what + 'は 0 より大きい数で入力してください（' + value + '）。');
    }
    return n;
  }

  function wholeBp(value) {
    var n = positive(value, '断片の長さ (bp) ');
    if (Math.abs(n - Math.round(n)) > 1e-9) {
      throw GibsonError('断片の長さ (bp) は整数で入力してください（' + value + '）。');
    }
    return Math.round(n);
  }

  function round(value, decimals) {
    var factor = Math.pow(10, decimals);
    return Math.round(value * factor) / factor;
  }

  /* Rounded to N significant figures. pmol in a Gibson reaction is often
   * 0.0x, where a fixed number of decimals either loses the value or prints
   * noise, so the number of digits follows the magnitude. */
  function significant(value, digits) {
    if (!value) { return 0; }
    var scale = digits - 1 - Math.floor(Math.log(Math.abs(value)) / Math.LN10);
    return round(value, scale);
  }

  // =====================================================================
  // ng / pmol / microlitres
  // =====================================================================

  function pmolFromNg(ng, bp) {
    return positive(ng, 'DNA の量 (ng) ') * 1000
      / (wholeBp(bp) * DS_DNA_G_PER_MOL_PER_BP);
  }

  function ngFromPmol(pmol, bp) {
    return positive(pmol, 'DNA の量 (pmol) ')
      * wholeBp(bp) * DS_DNA_G_PER_MOL_PER_BP / 1000;
  }

  function volumeUl(ng, concentration) {
    return positive(ng, 'DNA の量 (ng) ')
      / positive(concentration, '濃度 (ng/µL) ');
  }

  /* How much of each fragment to take.
   *
   * basePmol is the amount that molar ratio 1 stands for. Which fragment is
   * the "1" is the user's choice -- usually the vector -- and this function
   * does not care. No default is supplied for either basePmol or the ratios:
   * that is the point (see the spec, section 2).
   *
   * concentration may be left out, in which case that row gets no volume and
   * neither does the total. A partial answer is still worth showing. */
  function plan(spec) {
    var settings = spec || {};
    var base = positive(settings.basePmol, 'モル比 1 にあたる量 (pmol) ');
    var fragments = settings.fragments || [];
    if (!fragments.length) {
      throw GibsonError('断片が 1 つも入力されていません。');
    }

    var notes = [];
    var rawPmol = 0, rawNg = 0, roundedUl = 0, everyRowHasVolume = true;

    var rows = fragments.map(function (fragment, index) {
      var label = String(fragment.name || ('断片 ' + (index + 1)));
      var bp = wholeBp(fragment.bp);
      var ratio = positive(fragment.ratio, label + ' のモル比 ');
      var pmol = base * ratio;
      var ng = pmol * bp * DS_DNA_G_PER_MOL_PER_BP / 1000;
      rawPmol += pmol;
      rawNg += ng;

      var ul = null;
      if (fragment.concentration === undefined
          || fragment.concentration === null
          || fragment.concentration === '') {
        everyRowHasVolume = false;
      } else {
        ul = round(ng / positive(fragment.concentration,
                                 label + ' の濃度 (ng/µL) '), UL_DECIMALS);
        roundedUl += ul;
        if (ul < MIN_PIPETTE_UL) {
          notes.push(label + ' は ' + ul + ' µL です。' + MIN_PIPETTE_UL
            + ' µL 未満はピペットの誤差が大きくなります。'
            + '希釈するか、反応を大きくしてください。');
        }
      }

      return { name: label, bp: bp, ratio: ratio,
               pmol: significant(pmol, PMOL_DIGITS),
               ng: significant(ng, NG_DIGITS), ul: ul };
    });

    // pmol and ng are summed before rounding, because that is what is in the
    // tube chemically. Volumes are summed after rounding, because the rounded
    // number is what actually gets pipetted.
    var totalUl = everyRowHasVolume ? round(roundedUl, UL_DECIMALS) : null;

    var limit = settings.dnaVolumeLimitUl;
    if (totalUl !== null && limit !== undefined && limit !== null
        && limit !== '') {
      var allowed = positive(limit, 'DNA に使える体積 (µL) ');
      if (totalUl > allowed) {
        notes.push('DNA の合計が ' + totalUl + ' µL で、入れられる ' + allowed
          + ' µL を超えています。濃い液を用意するか、pmol を下げてください。');
      }
    }

    return { rows: rows,
             totals: { pmol: significant(rawPmol, PMOL_DIGITS),
                       ng: significant(rawNg, NG_DIGITS), ul: totalUl },
             notes: notes };
  }

  /* The reaction as it is laid out on the bench sheet: the vector first,
   * then the inserts, each with a hand-entered volume; the master mix and
   * water fill the rest of the reaction.
   *
   * Only concentration, length and the volumes are inputs; the total is
   * always TOTAL_UL.
   * Everything else is derived, and the guide above is checked, not applied.
   * pmol is rounded to 4 decimals here to match the bench sheet (0.0171,
   * 0.3996), where plan() uses 3 significant figures. */
  /* The name a fragment gets when the user gives it none: the vector, then
   * "Insert" for a single insert, or Insert A, B, C... (requested
   * 2026-10-11; numbers were used before). */
  function fragmentLabel(index, count) {
    if (index === 0) { return 'Vector'; }
    return count === 2 ? 'Insert' : 'Insert ' + String.fromCharCode(64 + index);
  }

  function reaction(spec) {
    var settings = spec || {};
    var fragments = settings.fragments || [];
    if (fragments.length < 2) {
      throw GibsonError('Vector と Insert を 1 本以上入力してください。');
    }
    var total = TOTAL_UL;
    var mixUl = round(total / MASTER_MIX_FOLD, UL_DECIMALS);

    var rawPmol = 0, dnaUl = 0;
    var rows = fragments.map(function (fragment, index) {
      var label = String(fragment.name || fragmentLabel(index, fragments.length));
      var concentration = positive(fragment.concentration, label + ' の濃度 (ng/µL) ');
      var bp = wholeBp(fragment.bp);
      var ul = positive(fragment.ul, label + ' の Volume (µL) ');
      var pmolPerUl = concentration * 1000 / (bp * DS_DNA_G_PER_MOL_PER_BP);
      rawPmol += ul * pmolPerUl;
      dnaUl += ul;
      return { name: label, concentration: concentration, bp: bp, ul: ul,
               pmolPerUl: round(pmolPerUl, 4), pmol: round(ul * pmolPerUl, 4),
               ng: significant(ul * concentration, NG_DIGITS),
               rawNg: ul * concentration };
    });

    dnaUl = round(dnaUl, UL_DECIMALS);
    // More DNA than the space left by the master mix is not a reaction that
    // can be set up, so it is refused rather than judged (user, 2026-10-10).
    var dnaRoom = round(total - mixUl, UL_DECIMALS);
    if (dnaUl > dnaRoom) {
      throw GibsonError('Insert と Vector の合計が ' + dnaUl + ' µL で、' + dnaRoom
        + ' µL（全量 ' + total + ' µL − G.A. M ' + mixUl + ' µL）を超えるため計算しません。'
        + 'Volume を減らすか、濃い DNA を用意してください。');
    }
    var waterUl = round(total - mixUl - dnaUl, UL_DECIMALS);
    var totalPmol = round(rawPmol, 4);
    var vector = rows[0];

    function within(v, range) { return v >= range[0] && v <= range[1]; }

    var checks = [{ label: 'Vector の量', target: VECTOR_NG[0] + '–' + VECTOR_NG[1] + ' ng',
                    value: vector.ng + ' ng', ok: within(vector.ng, VECTOR_NG) }];
    // Per-cell verdicts, so the screen can colour each number by the check
    // it belongs to (blue = meets the guide, red = does not; 2026-10-10).
    vector.ngOk = checks[0].ok;
    rows.slice(1).forEach(function (row) {
      var ratio = round(row.rawNg / vector.rawNg, 2);
      var short = row.bp < SHORT_INSERT_BP;
      row.ratio = ratio;
      checks.push({
        label: row.name + ' / Vector（ng 比、' + (short ? SHORT_INSERT_BP + ' bp 未満'
                                                   : SHORT_INSERT_BP + ' bp 以上') + '）',
        target: short ? '×' + INSERT_RATIO_SHORT_MIN + ' 以上'
                      : '×' + INSERT_RATIO_LONG[0] + '–' + INSERT_RATIO_LONG[1],
        value: '×' + ratio,
        ok: short ? ratio >= INSERT_RATIO_SHORT_MIN : within(ratio, INSERT_RATIO_LONG) });
      row.ratioOk = checks[checks.length - 1].ok;
    });
    var range = null;
    TOTAL_PMOL_RANGES.forEach(function (r) {
      if (rows.length >= r.min && rows.length <= r.max) { range = r.pmol; }
    });
    checks.push({ label: '合計 pmol（' + rows.length + ' fragment）',
                  target: range ? range[0] + '–' + range[1] + ' pmol' : '目安なし',
                  value: totalPmol + ' pmol',
                  ok: range ? within(totalPmol, range) : null });
    var pmolOk = range ? within(totalPmol, range) : null;
    rows.forEach(function (row) { row.ulOk = row.ul >= MIN_AMOUNT_UL; });
    var smallest = rows.reduce(function (a, b) { return b.ul < a.ul ? b : a; });
    checks.push({ label: 'ピペットで取れる量（最小の Volume）',
                  target: MIN_AMOUNT_UL + ' µL 以上',
                  value: smallest.ul + ' µL（' + smallest.name + '）',
                  ok: smallest.ul >= MIN_AMOUNT_UL });

    // No under-0.5 uL note here: the user found it unnecessary for this
    // table (2026-10-10). plan() still has it.
    var notes = [];
    rows.forEach(function (row) { delete row.rawNg; });

    return { rows: rows, mixUl: mixUl, waterUl: waterUl, dnaUl: dnaUl,
             totals: { pmol: totalPmol, ul: total, pmolOk: pmolOk },
             checks: checks, notes: notes };
  }

  /* The guide as a list, for showing on screen. Built from the same
   * constants reaction() checks against, so the two cannot drift apart. */
  // The DNA room, total volume and pipette minimum are still checked, but
  // the user asked for them not to be listed (2026-10-10).
  function guide() {
    return [
      { label: 'Vector の量', target: VECTOR_NG[0] + '–' + VECTOR_NG[1] + ' ng' },
      { label: 'Insert / Vector（ng 比）、Insert ' + SHORT_INSERT_BP + ' bp 以上',
        target: '×' + INSERT_RATIO_LONG[0] + '–' + INSERT_RATIO_LONG[1] },
      { label: 'Insert / Vector（ng 比）、Insert ' + SHORT_INSERT_BP + ' bp 未満',
        target: '×' + INSERT_RATIO_SHORT_MIN + ' 以上' }
    ].concat(TOTAL_PMOL_RANGES.map(function (r) {
      return { label: '合計 pmol、' + r.min + '–' + r.max + ' fragment（Vector + Insert）',
               target: r.pmol[0] + '–' + r.pmol[1] + ' pmol' };
    }));
  }

  /* The second way of filling the table (requested 2026-10-10): the user
   * enters the vector as usual -- Conc, Length and Volume -- and, for each
   * insert, how many times the vector's ng to put in. The vector's ng is
   * Volume x Conc (until the evening of 2026-10-10 it was typed in as ng;
   * the user asked for it to come from the table). Insert volumes follow
   * from their concentrations, rounded to the
   * 0.01 uL that is pipetted, and the result goes through reaction() so it
   * is judged exactly like a hand-entered table. Because of that rounding
   * the ng actually added differs a little from the target; both are
   * returned (targets) so the user can check the arithmetic. */
  function amountsFromRatio(spec) {
    var settings = spec || {};
    var fragments = settings.fragments || [];
    if (fragments.length < 2) {
      throw GibsonError('Vector と Insert を 1 本以上入力してください。');
    }
    var vectorUl = positive(fragments[0].ul, String(fragments[0].name || 'Vector')
                            + ' の Volume (µL) ');
    var vectorNg = vectorUl * positive(fragments[0].concentration,
                                       String(fragments[0].name || 'Vector') + ' の濃度 (ng/µL) ');
    var targets = fragments.map(function (fragment, index) {
      var label = String(fragment.name || fragmentLabel(index, fragments.length));
      var ratio = index === 0 ? 1 : positive(fragment.ratio, label + ' の倍数 ');
      var concentration = positive(fragment.concentration, label + ' の濃度 (ng/µL) ');
      var ng = vectorNg * ratio;
      // The vector's own volume is taken as typed, not recomputed from ng.
      var ul = index === 0 ? vectorUl : round(ng / concentration, UL_DECIMALS);
      if (ul <= 0) {
        throw GibsonError(label + ' は ' + significant(ng, NG_DIGITS) + ' ng で、0.01 µL '
          + '未満になります。DNA を希釈してください。');
      }
      return { name: label, ratio: ratio, targetNg: significant(ng, NG_DIGITS), ul: ul };
    });
    var result = reaction({ fragments: fragments.map(function (fragment, i) {
      return { name: fragment.name, concentration: fragment.concentration,
               bp: fragment.bp, ul: targets[i].ul };
    }) });
    result.targets = targets.map(function (target, i) {
      var row = result.rows[i];
      return { name: row.name, targetNg: target.targetNg, ul: row.ul, actualNg: row.ng,
               targetRatio: i === 0 ? null : target.ratio,
               actualRatio: i === 0 ? null : row.ratio };
    });
    result.vectorNg = significant(vectorNg, NG_DIGITS);
    return result;
  }

  /* Amounts that meet the guide, offered for the user to accept or not
   * (requested 2026-10-10). The user still decides: nothing is written into
   * the table until they press the button that copies these in.
   *
   * Every candidate is put through reaction() with its volumes rounded to
   * 0.01 uL -- what will actually be pipetted -- so a recommendation can
   * never pass here and then fail the checks on screen.
   *
   * Pass 1 looks for amounts meeting every check, over every volume that
   * can be pipetted (0.01 uL steps), so that if such amounts exist they are
   * found. Until 2026-10-10 evening it tried the vector in 1 ng steps and
   * moved all insert ratios together, which could miss a solution; the user
   * pointed this out. For each vector volume giving 50-100 ng, the volume
   * window of each insert that meets its ratio is worked out exactly; the
   * smallest volumes give the least pmol and the least water used, so if
   * they break the upper pmol limit or the reaction volume, no choice of
   * insert volumes can work for that vector volume. Otherwise inserts are
   * moved toward the middle of their ratio ranges as far as the limits
   * allow, or raised (most pmol per uL first) to reach a lower pmol limit.
   * Among the solutions it prefers the vector nearest 75 ng, then insert
   * ratios nearest the middle of x2-3 (for "x5 or more", nearest x5).
   * Every volume is at least MIN_AMOUNT_UL (0.22 uL), the smallest the
   * user's pipette takes; that is one of reaction()'s checks, so the search
   * starts its windows there.
   *
   * Short inserts often make pass 1 impossible: 50 ng of a 5 kb vector with
   * a 300 bp insert at x2 is already 0.514 pmol. The bench sheet then keeps
   * the vector amount and the total pmol and lets the ratio fall (x1.38 on
   * 2026-10-10), so pass 2 does the same: every insert ratio is scaled down
   * from the bottom of its range by one common factor, as little as
   * possible. The user confirmed this priority on 2026-10-10. */
  var CENTI = 100;               // volumes are searched in 0.01 uL units
  var RATIO_ROUND_SLACK = 0.005; // reaction() rounds the ratio to 2 decimals
  var T_STEPS = 50;              // positions tried between the window ends

  function recommend(spec) {
    var settings = spec || {};
    var fragments = settings.fragments || [];
    // Validate every input once, with the same messages reaction() gives.
    reaction({ fragments: fragments.map(function (f) {
      return { name: f.name, concentration: f.concentration, bp: f.bp, ul: 1 };
    }) });
    var total = TOTAL_UL;
    var conc = fragments.map(function (f) { return parseFloat(f.concentration); });
    var bps = fragments.map(function (f) { return wholeBp(f.bp); });
    var short = bps.map(function (bp) { return bp < SHORT_INSERT_BP; });
    // pmol contributed by 0.01 uL of each fragment
    var perCenti = conc.map(function (c, i) {
      return c * 1000 / (bps[i] * DS_DNA_G_PER_MOL_PER_BP) / CENTI;
    });
    var dnaLimit = Math.floor((total - round(total / MASTER_MIX_FOLD, UL_DECIMALS))
                              * CENTI + 1e-6);
    var pmolRange = null;
    TOTAL_PMOL_RANGES.forEach(function (r) {
      if (fragments.length >= r.min && fragments.length <= r.max) { pmolRange = r.pmol; }
    });

    // Over the DNA room reaction() refuses to calculate, so such volumes are
    // not candidates at all.
    function evaluate(centis) {
      if (volumeOf(centis) > dnaLimit) { return null; }
      return reaction({ fragments: fragments.map(function (f, i) {
        return { name: f.name, concentration: f.concentration, bp: f.bp,
                 ul: round(centis[i] / CENTI, UL_DECIMALS) };
      }) });
    }
    function attempt(vectorNg, ratioFor) {
      var centis = fragments.map(function (f, i) {
        var ng = i === 0 ? vectorNg : vectorNg * ratioFor(i);
        return Math.round(ng / conc[i] * CENTI);
      });
      if (centis.some(function (c) { return c <= 0; })) { return null; }
      return evaluate(centis);
    }
    function fails(result, skipRatios) {
      return result.checks.filter(function (c, i) {
        var isRatio = i > 0 && i < result.rows.length;
        return c.ok === false && !(skipRatios && isRatio);
      });
    }
    function pmolOf(centis) {
      return centis.reduce(function (sum, c, i) { return sum + c * perCenti[i]; }, 0);
    }
    function volumeOf(centis) {
      return centis.reduce(function (sum, c) { return sum + c; }, 0);
    }
    var middle = (VECTOR_NG[0] + VECTOR_NG[1]) / 2;
    var inserts = fragments.map(function (f, i) { return i; }).slice(1);

    // ---- pass 1: every check met, over all 0.01 uL volumes
    var best = null, bestScore = Infinity, fewest = null, overRoom = false;
    var maxVectorCenti = Math.ceil(VECTOR_NG[1] / conc[0] * CENTI) + 1;
    var minCenti = Math.round(MIN_AMOUNT_UL * CENTI);
    for (var cv = minCenti; cv <= Math.min(maxVectorCenti, dnaLimit); cv++) {
      var vectorNg = significant(cv * conc[0] / CENTI, NG_DIGITS);
      if (vectorNg < VECTOR_NG[0] || vectorNg > VECTOR_NG[1]) { continue; }
      var rawVectorNg = cv * conc[0] / CENTI;

      // The window of volumes meeting each insert's ratio.
      var lo = [], hi = [], aim = [], empty = false;
      inserts.forEach(function (i) {
        var minRatio = short[i] ? INSERT_RATIO_SHORT_MIN : INSERT_RATIO_LONG[0];
        var low = Math.max(minCenti, Math.ceil((minRatio - RATIO_ROUND_SLACK)
                                        * rawVectorNg / conc[i] * CENTI - 1e-9));
        var high = short[i] ? dnaLimit
          : Math.floor((INSERT_RATIO_LONG[1] + RATIO_ROUND_SLACK) * rawVectorNg
                       / conc[i] * CENTI - 1e-9);
        var target = short[i] ? low
          : Math.round((INSERT_RATIO_LONG[0] + INSERT_RATIO_LONG[1]) / 2
                       * rawVectorNg / conc[i] * CENTI);
        if (high < low) { empty = true; }
        lo[i] = low; hi[i] = high; aim[i] = Math.min(Math.max(target, low), high);
      });
      if (empty) { continue; }

      var base = fragments.map(function (f, i) { return i === 0 ? cv : lo[i]; });
      // Smallest insert volumes = least pmol and least volume.
      var baseResult = evaluate(base);
      if (!baseResult) { overRoom = true; continue; }
      var baseFailed = fails(baseResult, false);
      if (!fewest || baseFailed.length < fewest.failed.length) {
        fewest = { result: baseResult, failed: baseFailed };
      }
      if (pmolRange && round(pmolOf(base), 4) > pmolRange[1]) { continue; }

      // Walk from the smallest volumes toward the aim; keep the furthest
      // point that still fits. t = 1 is the aim itself.
      var chosen = null, chosenT = -1;
      for (var step = T_STEPS; step >= 0; step--) {
        var t = step / T_STEPS;
        var centis = base.map(function (c, i) {
          return i === 0 ? c : Math.round(lo[i] + (aim[i] - lo[i]) * t);
        });
        var pmol = round(pmolOf(centis), 4);
        if (volumeOf(centis) > dnaLimit) { continue; }
        if (pmolRange && (pmol > pmolRange[1] || pmol < pmolRange[0])) { continue; }
        chosen = centis; chosenT = t;
        break;
      }
      // A lower pmol limit not reached at the aim: raise the inserts that add
      // the most pmol per uL first, within their windows and the volume.
      if (!chosen && pmolRange) {
        var raised = base.slice();
        inserts.slice().sort(function (a, b) { return perCenti[b] - perCenti[a]; })
          .forEach(function (i) {
            while (round(pmolOf(raised), 4) < pmolRange[0] && raised[i] < hi[i]
                   && volumeOf(raised) < dnaLimit) { raised[i]++; }
          });
        var raisedPmol = round(pmolOf(raised), 4);
        if (raisedPmol >= pmolRange[0] && raisedPmol <= pmolRange[1]
            && volumeOf(raised) <= dnaLimit) { chosen = raised; chosenT = 0; }
      }
      if (!chosen) { continue; }

      var result = evaluate(chosen);
      if (!result || fails(result, false).length) { continue; }   // rounding edge
      var score = Math.abs(vectorNg - middle) / middle + (1 - chosenT);
      if (score < bestScore) { best = result; bestScore = score; }
    }
    if (best) {
      return { satisfied: true, relaxed: false, result: best,
               message: '条件設定をすべて満たす Volume です。' };
    }

    // No vector volume at all gives 50-100 ng: say which way it is off,
    // since no change to the inserts can fix that.
    var vectorAtMin = significant(minCenti * conc[0] / CENTI, NG_DIGITS);
    if (vectorAtMin > VECTOR_NG[1]) {
      return { satisfied: false, relaxed: false, result: null,
               message: 'Vector は最小の ' + MIN_AMOUNT_UL + ' µL でも ' + vectorAtMin
                 + ' ng で、' + VECTOR_NG[1] + ' ng を超えます。Vector を希釈してください。' };
    }
    var vectorAtLimit = significant(dnaLimit * conc[0] / CENTI, NG_DIGITS);
    if (vectorAtLimit < VECTOR_NG[0]) {
      return { satisfied: false, relaxed: false, result: null,
               message: 'Vector は DNA に使える ' + round(dnaLimit / CENTI, UL_DECIMALS)
                 + ' µL を全部使っても ' + vectorAtLimit + ' ng で、' + VECTOR_NG[0]
                 + ' ng に届きません。濃い Vector を用意してください。' };
    }

    // ---- pass 2: keep the vector, total pmol and the volumes; lower the ratios
    for (var k = 99; k >= 1; k--) {
      var factor = k / 100, pick = null, pickScore = Infinity;
      var lowered = function (i) {
        return (short[i] ? INSERT_RATIO_SHORT_MIN : INSERT_RATIO_LONG[0]) * factor;
      };
      for (var w = VECTOR_NG[0]; w <= VECTOR_NG[1]; w++) {
        var candidate = attempt(w, lowered);
        if (!candidate || fails(candidate, true).length) { continue; }
        var s2 = Math.abs(w - middle) / middle;
        if (s2 < pickScore) { pick = candidate; pickScore = s2; }
      }
      if (pick) {
        return { satisfied: false, relaxed: true, result: pick,
                 message: 'Insert / Vector の倍数を目安どおりにすると、ほかの条件'
                   + '（合計 pmol・DNA の合計 ' + round(dnaLimit / CENTI, UL_DECIMALS)
                   + ' µL 以下・ピペットで取れる量）を外れます。'
                   + 'それらと Vector の量を守り、'
                   + '倍数を目安の下限から最小限だけ下げた Volume です。' };
      }
    }

    return { satisfied: false, relaxed: false, result: null,
             message: '条件を満たす Volume は見つかりませんでした。満たせない条件: '
               + (fewest ? fewest.failed.map(function (c) { return c.label; }).join('、')
                  : overRoom ? 'Insert と Vector の合計を '
                    + round(dnaLimit / CENTI, UL_DECIMALS) + ' µL 以下にすること'
                  : '（候補を作れませんでした）') + '。' };
  }

  // =====================================================================
  // Tm
  // =====================================================================

  function clean(text) {
    return String(text === undefined || text === null ? '' : text)
      .toUpperCase().replace(/U/g, 'T').replace(/[^A-Z]/g, '');
  }

  /* Tm, GC and length of one oligo.
   *
   * The Takara formula is stated for 25 mer or shorter, and a Gibson primer
   * with its overlap is routinely 40-60 mer. The value is still returned --
   * hiding it would only send the user to a different calculator -- with
   * formulaApplies false and a note saying so. Switching to a second formula
   * past 25 mer was considered and rejected: two formulas means never knowing
   * which number is on screen. */
  function tmOf(sequence) {
    var seq = clean(sequence);
    if (!seq.length) { throw GibsonError('配列を入力してください。'); }
    var pb = engine();
    var tm = pb.tmTakara(seq);          // throws on anything but A/C/G/T
    var applies = pb.tmFormulaApplies(seq.length);
    return { sequence: seq, length: seq.length, tm: tm,
             gcPercent: round(pb.gcPercent(seq), GC_DECIMALS),
             formulaApplies: applies,
             note: applies ? null
               : seq.length + ' mer は簡易式の適用範囲（'
                 + pb.TM_FORMULA_MAX_LEN + ' mer 以下）を超えています。'
                 + '目安として見てください。' };
  }

  // One sequence per line. Spaces and numbers inside a line are stripped, so
  // a block copied out of ApE pastes in as it is.
  function tmRows(text) {
    return String(text === undefined || text === null ? '' : text)
      .split('\n')
      .map(function (line) { return clean(line); })
      .filter(function (seq) { return seq.length > 0; })
      .map(tmOf);
  }

  function tmDifference(a, b) {
    return Math.abs(tmOf(a).tm - tmOf(b).tm);
  }

  /* PCR conditions for the primers in the Tm block, from the enzyme's
   * protocol the user pasted on 2026-10-10 (Takara Bio product page
   * unitid=U100004616). Only what that text states is returned: there is no
   * extension time or cycle count in it, so none is made up here.
   *
   * The annealing time rule is PB.annealingTimeSec -- the same rule Primer
   * Builder applies -- and, as there, the reaction takes the longest time any
   * one primer needs (PB.annealingTimeForPair).
   *
   * The 70 C two-step hint is judged only for primers within the formula's
   * 25 mer range: a 50 mer Gibson primer always "scores" over 70 by the
   * formula, which says nothing about the primer. */
  var TWO_STEP_TM_C = 70;
  var DENATURE_98_SEC = '5–10 sec', DENATURE_94_SEC = '10–15 sec';

  function pcrConditions(rows) {
    if (!rows || !rows.length) { throw GibsonError('配列を入力してください。'); }
    var pb = engine();
    var perPrimer = rows.map(function (row) {
      return pb.annealingTimeSec(row.length, row.tm);
    });
    var hot = [];
    rows.forEach(function (row, i) {
      if (row.formulaApplies && row.tm >= TWO_STEP_TM_C) { hot.push(i + 1); }
    });
    var notes = ['3 step でスメアになる場合は、2 step での反応を試す。'];
    if (hot.length) {
      notes.unshift('#' + hot.join('、#') + ' の Tm が ' + TWO_STEP_TM_C
        + '℃ 以上です。2 step での反応を試してください。');
    }
    return { mode: '3 step',
             denature: { tempC: 98, time: DENATURE_98_SEC,
                         alternative: '94℃ の場合は ' + DENATURE_94_SEC },
             anneal: { tempC: pb.ANNEAL_TEMP_C,
                       sec: Math.max.apply(null, perPrimer) },
             perPrimerSec: perPrimer, twoStepRows: hot, notes: notes };
  }

  return {
    DS_DNA_G_PER_MOL_PER_BP: DS_DNA_G_PER_MOL_PER_BP,
    normalizeNumberText: normalizeNumberText,
    MIN_PIPETTE_UL: MIN_PIPETTE_UL,
    MASTER_MIX_FOLD: MASTER_MIX_FOLD, TOTAL_UL: TOTAL_UL, MIN_AMOUNT_UL: MIN_AMOUNT_UL,
    round: round, significant: significant,
    pmolFromNg: pmolFromNg, ngFromPmol: ngFromPmol, volumeUl: volumeUl,
    plan: plan, reaction: reaction, recommend: recommend, guide: guide,
    fragmentLabel: fragmentLabel,
    amountsFromRatio: amountsFromRatio,
    tmOf: tmOf, tmRows: tmRows, tmDifference: tmDifference,
    pcrConditions: pcrConditions
  };
})();
if (typeof module !== 'undefined' && module.exports) { module.exports = GIBSON; }
