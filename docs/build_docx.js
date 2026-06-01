// Генератор docs/algorithms.docx по ГОСТ 2.105 / 7.32 / РД 50-34.698-90.
// Запуск:  node docs/build_docx.js
// Требует:  npm install -g docx  (или локально)

const fs   = require('fs');
const path = require('path');

// Подхватываем глобально установленный пакет.
const GLOBAL_NM = (process.env.NPM_GLOBAL_NM
                || 'C:/Users/Sasha Bray/AppData/Roaming/npm/node_modules');
const docx = require(path.join(GLOBAL_NM, 'docx'));

const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, Header, Footer,
  AlignmentType, LevelFormat, HeadingLevel, BorderStyle, WidthType, ShadingType,
  PageNumber, PageBreak, TabStopType, TabStopPosition,
  Math: MathBlock, MathRun, MathFraction, MathRadical, MathSuperScript, MathSubScript,
  MathSubSuperScript, MathRoundBrackets, MathSquareBrackets, MathAngledBrackets, MathSum
} = docx;

// ---------------------------------------------------------------------------
// Геометрия страницы A4 + поля по ГОСТ 7.32 (30 / 15 / 20 / 20 мм)
// 1440 DXA = 1 дюйм = 25.4 мм
// ---------------------------------------------------------------------------
const mm = n => Math.round(n * 1440 / 25.4);

const PAGE = {
  width:  mm(210),
  height: mm(297),
  marginLeft:   mm(30),
  marginRight:  mm(15),
  marginTop:    mm(20),
  marginBottom: mm(20),
};
const CONTENT_W = PAGE.width - PAGE.marginLeft - PAGE.marginRight;
const FLI = mm(12.5);             // первая строка — 12.5 мм
const TIMES = 'Times New Roman';
const MONO  = 'Consolas';

// Размер — в half-points: 28 = 14 pt
const SIZE_BODY = 28;
const SIZE_CODE = 22;             // 11 pt
const SIZE_CAPT = 24;             // 12 pt

// ---------------------------------------------------------------------------
// Математика: helper'ы поверх docx-js Math API.
// mathPara(...nodes) — один центрированный абзац с одной формулой.
// Каждый узел может быть либо строкой (превратится в MathRun),
// либо объектом Math* (используется как есть).
// ---------------------------------------------------------------------------
function _m(x) { return typeof x === 'string' ? new MathRun(x) : x; }
function _arr(x) { return (Array.isArray(x) ? x : [x]).map(_m); }
function mathPara(...nodes) {
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { line: 360, lineRule: 'auto', before: 160, after: 160 },
    indent: { firstLine: 0 },
    children: [ new MathBlock({ children: nodes.map(_m) }) ],
  });
}
// Несколько строк в одной формуле — каждый row рендерится отдельным абзацем.
function mathParas(rows) {
  return rows.map(r => mathPara(...(Array.isArray(r) ? r : [r])));
}
function frac(num, den) {
  return new MathFraction({ numerator: _arr(num), denominator: _arr(den) });
}
function sup(base, exp) {
  return new MathSuperScript({ children: _arr(base), superScript: _arr(exp) });
}
function sub(base, idx) {
  return new MathSubScript({ children: _arr(base), subScript: _arr(idx) });
}
function subSup(base, idx, exp) {
  return new MathSubSuperScript({
    children: _arr(base), subScript: _arr(idx), superScript: _arr(exp),
  });
}
function rad(content) {
  return new MathRadical({ children: _arr(content) });
}
function paren(content) {
  return new MathRoundBrackets({ children: _arr(content) });
}
function bra(content) {
  return new MathSquareBrackets({ children: _arr(content) });
}

// ---------------------------------------------------------------------------
// Стилевые помощники
// ---------------------------------------------------------------------------
function run(text, opts = {}) {
  return new TextRun({ text, font: TIMES, size: SIZE_BODY, ...opts });
}
function body(text, opts = {}) {
  const runs = Array.isArray(text) ? text : [run(text)];
  return new Paragraph({
    spacing: { line: 360, lineRule: 'auto', after: 60 },
    indent: { firstLine: FLI },
    alignment: AlignmentType.JUSTIFIED,
    children: runs,
    ...opts,
  });
}
function bodyNoIndent(text, opts = {}) {
  return body(text, { indent: { firstLine: 0 }, ...opts });
}
// Заголовок раздела (1, 2, ...) — прописными, по центру, жирный, по ГОСТ 7.32
function h1(num, text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    alignment: AlignmentType.CENTER,
    spacing: { before: 480, after: 240, line: 360, lineRule: 'auto' },
    pageBreakBefore: true,
    children: [
      new TextRun({ text: `${num} ${text}`.toUpperCase(),
                    font: TIMES, size: SIZE_BODY, bold: true }),
    ],
  });
}
// Подзаголовок 1.1, 1.2 — с красной строки, жирный, без точки
function h2(num, text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_2,
    indent: { firstLine: FLI },
    spacing: { before: 280, after: 140, line: 360, lineRule: 'auto' },
    children: [
      new TextRun({ text: `${num} ${text}`, font: TIMES, size: SIZE_BODY, bold: true }),
    ],
  });
}
// Заголовок без номера (Содержание, Введение, Приложение А) — ПРОПИСНЫМИ, по центру
function h1Unnumbered(text, opts = {}) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    alignment: AlignmentType.CENTER,
    spacing: { before: 480, after: 360, line: 360, lineRule: 'auto' },
    pageBreakBefore: opts.pageBreakBefore !== false,
    children: [
      new TextRun({ text: text.toUpperCase(), font: TIMES, size: SIZE_BODY, bold: true }),
    ],
  });
}
// Многострочный псевдокод (моноширинный, без отступа красной строки)
function code(text) {
  const lines = text.split('\n');
  const children = [];
  lines.forEach((line, i) => {
    if (i > 0) children.push(new TextRun({ break: 1 }));
    children.push(new TextRun({ text: line.replace(/^\s+/, s => ' '.repeat(s.length)),
                                font: MONO, size: SIZE_CODE }));
  });
  return new Paragraph({
    spacing: { line: 276, lineRule: 'auto', before: 120, after: 120 },
    indent: { left: mm(7), firstLine: 0 },
    alignment: AlignmentType.LEFT,
    children,
  });
}
// Подпись таблицы — «Таблица N — Название», слева, без красной строки, перед таблицей
function tableCaption(num, title) {
  return new Paragraph({
    spacing: { before: 200, after: 100, line: 360, lineRule: 'auto' },
    alignment: AlignmentType.LEFT,
    indent: { firstLine: 0 },
    children: [
      new TextRun({ text: `Таблица ${num} — ${title}`, font: TIMES, size: SIZE_BODY }),
    ],
  });
}
// Маркированный пункт (тире — традиционный маркер ГОСТ)
function dash(text) {
  return new Paragraph({
    spacing: { line: 360, lineRule: 'auto', after: 60 },
    indent: { left: mm(7.5), hanging: mm(5) },
    alignment: AlignmentType.JUSTIFIED,
    children: [ new TextRun({ text: '— ' + text, font: TIMES, size: SIZE_BODY }) ],
  });
}

// ---------------------------------------------------------------------------
// Таблицы — двухколоночные «Параметр / Описание», тонкая граница, паддинг
// ---------------------------------------------------------------------------
function thinBorder() {
  const b = { style: BorderStyle.SINGLE, size: 4, color: '000000' };
  return { top: b, bottom: b, left: b, right: b };
}
function tcell(text, opts = {}) {
  const w = opts.width;
  const isHeader = !!opts.header;
  return new TableCell({
    borders: thinBorder(),
    width: { size: w, type: WidthType.DXA },
    margins: { top: 80, bottom: 80, left: 120, right: 120 },
    shading: isHeader ? { fill: 'E8E8E8', type: ShadingType.CLEAR } : undefined,
    verticalAlign: 'top',
    children: [ new Paragraph({
      spacing: { line: 276, lineRule: 'auto' },
      alignment: AlignmentType.LEFT,
      indent: { firstLine: 0 },
      children: [ new TextRun({ text, font: TIMES, size: SIZE_CAPT, bold: isHeader }) ],
    }) ],
  });
}
function table2col(headers, rows, weights = [0.35, 0.65]) {
  const w1 = Math.round(CONTENT_W * weights[0]);
  const w2 = CONTENT_W - w1;
  const headRow = new TableRow({
    tableHeader: true,
    children: [ tcell(headers[0], { width: w1, header: true }),
                tcell(headers[1], { width: w2, header: true }) ],
  });
  const dataRows = rows.map(r => new TableRow({
    children: [ tcell(r[0], { width: w1 }), tcell(r[1], { width: w2 }) ],
  }));
  return new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: [w1, w2],
    rows: [headRow, ...dataRows],
  });
}
function table3col(headers, rows, weights = [0.22, 0.30, 0.48]) {
  const w1 = Math.round(CONTENT_W * weights[0]);
  const w2 = Math.round(CONTENT_W * weights[1]);
  const w3 = CONTENT_W - w1 - w2;
  const headRow = new TableRow({
    tableHeader: true,
    children: [ tcell(headers[0], { width: w1, header: true }),
                tcell(headers[1], { width: w2, header: true }),
                tcell(headers[2], { width: w3, header: true }) ],
  });
  const dataRows = rows.map(r => new TableRow({
    children: [ tcell(r[0], { width: w1 }),
                tcell(r[1], { width: w2 }),
                tcell(r[2], { width: w3 }) ],
  }));
  return new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: [w1, w2, w3],
    rows: [headRow, ...dataRows],
  });
}
function table4col(headers, rows, weights = [0.18, 0.30, 0.12, 0.40]) {
  const w1 = Math.round(CONTENT_W * weights[0]);
  const w2 = Math.round(CONTENT_W * weights[1]);
  const w3 = Math.round(CONTENT_W * weights[2]);
  const w4 = CONTENT_W - w1 - w2 - w3;
  const cells = headers.map((h, i) =>
    tcell(h, { width: [w1,w2,w3,w4][i], header: true }));
  const headRow = new TableRow({ tableHeader: true, children: cells });
  const dataRows = rows.map(r => new TableRow({
    children: r.map((c, i) => tcell(c, { width: [w1,w2,w3,w4][i] })),
  }));
  return new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: [w1, w2, w3, w4],
    rows: [headRow, ...dataRows],
  });
}

// ---------------------------------------------------------------------------
// Титульный лист
// ---------------------------------------------------------------------------
function blank() {
  return new Paragraph({ spacing: { before: 0, after: 0 },
                         children: [ new TextRun('') ] });
}
function titleLine(text, size = 28, bold = true, opts = {}) {
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: opts.before ?? 0, after: opts.after ?? 0,
               line: 360, lineRule: 'auto' },
    indent: { firstLine: 0 },
    children: [ new TextRun({ text, font: TIMES, size, bold }) ],
  });
}

const titlePage = [
  titleLine('Synchrotron Particle Tracker (SynchrotronSim)',
            SIZE_BODY, true, { before: 600 }),
  blank(),
  titleLine('РАСЧЁТНО-ПРОГРАММНЫЙ КОМПЛЕКС', SIZE_BODY, false),
  titleLine('моделирования траекторий частиц в магнитной структуре синхротрона',
            SIZE_BODY, false, { after: 1200 }),
  blank(), blank(), blank(),
  titleLine('ОПИСАНИЕ АЛГОРИТМОВ', 36, true, { before: 600, after: 240 }),
  blank(),
  titleLine('Оформлено по РД 50-34.698-90', 26, false),
  titleLine('(методические указания к ГОСТ 34, раздел «Описание алгоритма»)',
            26, false, { after: 1800 }),
  blank(), blank(), blank(), blank(), blank(), blank(),
  titleLine('2026', SIZE_BODY, false, { before: 1800 }),
];

// ---------------------------------------------------------------------------
// Содержание
// ---------------------------------------------------------------------------
// Простое — без авто-TOC: сохраняемое и при редактировании в Word.
function tocItem(num, text, page) {
  const tabs = [{ type: TabStopType.RIGHT, position: CONTENT_W - mm(2),
                  leader: 'dot' }];
  return new Paragraph({
    spacing: { line: 360, lineRule: 'auto', after: 60 },
    indent: { firstLine: 0 },
    tabStops: tabs,
    children: [
      new TextRun({ text: `${num} ${text}`, font: TIMES, size: SIZE_BODY }),
      new TextRun({ text: '\t' + page, font: TIMES, size: SIZE_BODY }),
    ],
  });
}

// Состав документа (для содержания)
const TOC_ENTRIES = [
  ['',  'Введение'],
  ['1', 'Интегратор Бориса — релятивистский симплектический интегратор'],
  ['2', 'Трилинейная интерполяция поля с учётом YZ-симметрии'],
  ['3', 'Двумерный bucket-индекс магнитов'],
  ['4', 'Вклад магнита с прескринингом по охватывающей сфере'],
  ['5', 'Детектор пересечения финишной плоскости'],
  ['6', 'Расстояние от точки до полилинии проектной орбиты'],
  ['7', 'Составная функция потерь оптимизатора'],
  ['8', 'Стохастический оптимизатор (hill-climber c киками)'],
  ['9', 'Генерация решётки из LTE-описания'],
  ['10','Сборка таблицы конфигурации и суб-решёток'],
  ['11','Анализ огибающей и бетатронных частот по FFT'],
  ['А', 'Приложение А. Файлы конфигурации'],
  ['Б', 'Приложение Б. Файлы данных'],
];

// ---------------------------------------------------------------------------
// Введение
// ---------------------------------------------------------------------------
const intro = [
  h1Unnumbered('Введение', { pageBreakBefore: true }),
  body('Настоящий документ содержит описание вычислительных моделей, реализованных в составе программ SynchrotronTracker.exe и SynchrotronOptimizer.exe. Документ оформлен в соответствии с методическими указаниями РД 50-34.698-90, разделом «Описание алгоритма», где для каждого алгоритма приводятся пять подразделов: назначение и характеристика, используемая информация, результаты решения, математическое описание, алгоритм решения.'),
  body('Источниками исходных данных служат файлы решётки в формате Elegant LTE (workspace/lte/), бинарные карты поля Ansys (workspace/fields_bin/), JSON-конфигурации (workspace/configs/) и CSV-таблицы базы данных (workspace/tables/). Результатами работы программ являются бинарные файлы траекторий (workspace/results/trajectories_exp{id}.bin) и журналы оптимизации в формате JSON.'),
  body('Все размерные величины приведены в системе СИ. Ссылки на исходный код имеют вид «файл:строка», например, src/Integrator.h:53. Документ описывает 11 алгоритмов и содержит два приложения со справочной информацией о конфигурационных файлах и форматах данных.'),
];

// ---------------------------------------------------------------------------
// 1. Boris pusher
// ---------------------------------------------------------------------------
const alg1 = [
  h1('1', 'Интегратор Бориса — релятивистский симплектический интегратор'),
  body('Исходный код: src/Integrator.h:53, метод Integrator::borisStep.'),

  h2('1.1', 'Назначение и характеристика'),
  body('Одношаговый релятивистский интегратор уравнений движения заряженной частицы в произвольном электромагнитном поле. Реализован по схеме Бориса с расщеплением drift–kick–drift. Обладает следующими свойствами:'),
  dash('симплектичен — сохраняет фазовый объём;'),
  dash('точно сохраняет модуль импульса |p|² для чисто магнитного поля (шаг поворота является ортогональным преобразованием по построению);'),
  dash('симметричен по времени, формальный порядок аппроксимации — O(dt²);'),
  dash('требует одного вычисления поля за шаг (против четырёх в методе Рунге–Кутты 4-го порядка), что даёт ускорение приблизительно в 4 раза по сравнению с РК4 при идентичном бюджете шагов.'),
  body('Литература: Boris J. P., «Relativistic plasma simulation», Proc. Fourth Conf. on Numerical Simulation of Plasmas, 1970; Birdsall C. K., Langdon A. B., «Plasma Physics via Computer Simulation», 1985, §4-3.'),

  h2('1.2', 'Используемая информация'),
  tableCaption('1.1', 'Входные данные алгоритма'),
  table4col(
    ['Обозначение', 'Смысл', 'Ед.', 'Источник'],
    [
      ['s.pos', 'координаты частицы на начало шага', 'м', 'ParticleState'],
      ['s.mom', 'импульс p = γ·m·v', 'кг·м/с', 'ParticleState'],
      ['m0',    'масса покоя', 'кг', 'particle_groups.csv'],
      ['q',     'заряд', 'Кл', 'particle_groups.csv'],
      ['lattice', 'структура с методом getTotalFields(pos)', '—', 'lattice_configs.csv'],
      ['dt',    'шаг интегрирования', 'с', 'experiments.csv'],
    ]),

  h2('1.3', 'Результаты решения'),
  tableCaption('1.2', 'Выходные данные алгоритма'),
  table2col(
    ['Обозначение', 'Смысл и место сохранения'],
    [
      ['s.pos', 'координаты частицы после шага (обновление по месту)'],
      ['s.mom', 'импульс частицы после шага (обновление по месту)'],
    ]),

  h2('1.4', 'Математическое описание'),
  body('Уравнения движения релятивистской заряженной частицы (уравнения Лоренца):'),
  mathPara(frac('dx', 'dt'), ' = v,    ', frac('dp', 'dt'), ' = q (E + v × B)'),
  body('Обозначения:'),
  dash('x(t) — координаты частицы в мировой системе, м.'),
  dash('v(t) — скорость частицы в мировой системе, м/с.'),
  dash('p(t) — релятивистский импульс, кг·м/с (величина, реально интегрируемая кодом; совпадает с классическим m·v только при малых скоростях).'),
  dash('q — заряд частицы, Кл; для электронов q = −1,602·10⁻¹⁹ Кл.'),
  dash('m — масса покоя, кг; для электронов m = 9,109·10⁻³¹ кг.'),
  dash('E, B — векторы электрического и магнитного полей в точке частицы, В/м и Тл.'),
  dash('t — физическое (лабораторное) время, с.'),
  body('Векторное произведение v × B есть магнитная сила Лоренца на единицу заряда: магнитное поле искривляет траекторию, но не совершает работу, и энергия сохраняется точно.'),
  body('Релятивистский импульс и Лоренц-фактор:'),
  mathPara('p = γ m v'),
  mathPara(
    'γ = ',
    frac('1', rad(['1 − ', frac(sup('|v|', '2'), sup('c', '2'))])),
    ' = ',
    rad(['1 + ', frac(sup('|p|', '2'), [sup('m', '2'), ' ', sup('c', '2')])]),
  ),
  body('Обозначения:'),
  dash('γ ≥ 1 — Лоренц-фактор, безразмерная величина, мера релятивистичности движения (γ = 1 в покое, γ → ∞ при |v| → c).'),
  dash('c = 299 792 458 м/с — скорость света (точное значение по определению метра).'),
  dash('|v|² = vx² + vy² + vz² — квадрат модуля скорости.'),
  dash('|p|² = px² + py² + pz² — квадрат модуля импульса.'),
  body('Вторая форма (через |p|) и используется в коде: γ восстанавливается прямо из p, не требуя хранения скорости. В текущей решётке все устройства магнитные, поэтому E ≡ 0; шаги E-толчка ниже фактически no-op, общая форма сохранена на случай появления электростатических элементов.'),
  body('Расщепление Бориса. Один шаг Δt состоит из семи подопераций над состоянием (xₙ, pₙ) → (xₙ₊₁, pₙ₊₁).'),
  body('Шаг 1 — полудрейф (координаты сдвигаются на половину шага по текущему импульсу):'),
  mathPara(sub('x', 'n+½'), ' = ', sub('x', 'n'), ' + ',
           frac(sub('p', 'n'), [sub('γ', 'n'), ' m']),
           ' · ', frac('Δt', '2')),
  dash('Δt — шаг интегрирования, с (берётся из experiments.csv.dt).'),
  dash('γ_n = √(1 + |p_n|² / (m²·c²)) — Лоренц-фактор на шаге n.'),
  body('Шаг 2 — опрос поля в середине шага:'),
  mathPara('(B, E) ← lattice.getTotalFields(', sub('x', 'n+½'), ')'),
  body('Это единственный вызов вычисления поля за шаг. РК4 потребовал бы четырёх.'),
  body('Шаг 3 — первый полу-E-толчок:'),
  mathPara(sup('p', '−'), ' = ', sub('p', 'n'), ' + ', frac('q E Δt', '2')),
  dash('p⁻ — промежуточный импульс после первого E-толчка. При E = 0 имеем p⁻ = p_n.'),
  body('Шаг 4 — пересчёт γ по новому импульсу:'),
  mathPara(sup('γ', '−'), ' = ', rad(['1 + ', frac(sup('|', sup('p', '−')), [sup('m', '2'), ' ', sup('c', '2')])])),
  // Замечание: формула выше — γ⁻ = √(1 + |p⁻|² / (m²c²)).
  body('Эта величина γ⁻ — «мгновенный» релятивистский фактор для последующего магнитного поворота.'),
  body('Шаг 5 — магнитный поворот. Сначала формируются вспомогательные векторы:'),
  mathPara('t = ', frac('q B Δt', ['2 ', sup('γ', '−'), ' m'])),
  mathPara('s = ', frac('2 t', ['1 + ', sup('|t|', '2')])),
  dash('t — половинный вектор поворота; его модуль равен tan(θ/2), где θ = qB·Δt/(γ⁻·m) — угол ларморовского поворота.'),
  dash('s — перенормированный вектор поворота, выбранный так, чтобы двухшаговое обновление ниже точно сохраняло модуль импульса.'),
  body('Затем импульс поворачивается в два этапа:'),
  mathPara('p′ = ', sup('p', '−'), ' + ', sup('p', '−'), ' × t'),
  mathPara(sup('p', '+'), ' = ', sup('p', '−'), ' + p′ × s,    |', sup('p', '+'), '| = |', sup('p', '−'), '| точно'),
  dash('p\' — промежуточный импульс после первого полуповорота.'),
  dash('p⁺ — импульс после полного магнитного поворота; равенство |p⁺| = |p⁻| справедливо как алгебраическое тождество — именно за счёт него Boris-pusher точно сохраняет энергию в чисто магнитном поле.'),
  body('Шаг 6 — второй полу-E-толчок:'),
  mathPara(sub('p', 'n+1'), ' = ', sup('p', '+'), ' + ', frac('q E Δt', '2')),
  body('Вместе с шагом 3 даёт центрированное E-ускорение: каждый полутолчок добавляет половину импульса q·E·Δt.'),
  body('Шаг 7 — второй полудрейф:'),
  mathPara(sub('x', 'n+1'), ' = ', sub('x', 'n+½'), ' + ',
           frac(sub('p', 'n+1'), [sub('γ', 'n+1'), ' m']),
           ' · ', frac('Δt', '2')),
  dash('γ_{n+1} = √(1 + |p_{n+1}|² / (m²·c²)) — Лоренц-фактор на новом шаге.'),
  body('После семи подопераций состояние продвинуто ровно на один Δt с порядком аппроксимации O(Δt²).'),

  h2('1.5', 'Алгоритм решения'),
  code(
    'Вход:  s = (pos, mom), m0, q, lattice, dt\n' +
    'Выход: s обновлено до шага n+1\n\n' +
    ' 1.  γ ← sqrt(1 + |mom|² / (m0²·c²))\n' +
    ' 2.  pos ← pos + mom · (dt/2) / (m0·γ)\n' +
    ' 3.  (B, E) ← lattice.getTotalFields(pos)\n' +
    ' 4.  p_minus ← mom + E · q · dt/2\n' +
    ' 5.  γ_minus ← sqrt(1 + |p_minus|² / (m0²·c²))\n' +
    ' 6.  t ← B · q · dt / (2·m0·γ_minus)\n' +
    ' 7.  s ← 2·t / (1 + |t|²)\n' +
    ' 8.  p_prime ← p_minus + cross(p_minus, t)\n' +
    ' 9.  p_plus  ← p_minus + cross(p_prime, s)\n' +
    '10.  mom ← p_plus + E · q · dt/2\n' +
    '11.  γ ← sqrt(1 + |mom|² / (m0²·c²))\n' +
    '12.  pos ← pos + mom · (dt/2) / (m0·γ)'),
  body('В ранних версиях существовали вспомогательные методы Integrator::driftStep и AcceleratorConfig::isDriftRegion для замены полного шага Бориса аналитическим дрейфом x += v·dt на участках без магнитов. Оба удалены в рефакторинге 2026-06: прескрининг по bucket-индексу (см. раздел 3) уже сводит шаг в пустой области к одному обращению к бакету и пустому циклу по кандидатам, и отдельная ветка для дрейфа практической пользы не давала.'),
  body('В текущей конфигурации все устройства магнитные, поэтому E ≡ 0 и шаги 3, 6 (полу-E-толчки) фактически no-op. Алгоритм оставлен в общем виде для совместимости при возможном появлении электрических элементов.'),
];

// ---------------------------------------------------------------------------
// 2. Trilinear field interpolation
// ---------------------------------------------------------------------------
const alg2 = [
  h1('2', 'Трилинейная интерполяция поля с учётом YZ-симметрии'),
  body('Исходный код: src/FieldMap.h:41, метод FieldMap::getInterpolatedField.'),

  h2('2.1', 'Назначение и характеристика'),
  body('Вычисление вектора магнитного (или электрического) поля в произвольной точке локальной системы координат магнита по заранее посчитанной на регулярной декартовой сетке карте поля (выход Ansys). Трилинейная интерполяция обеспечивает непрерывность C⁰ и совпадение со значением в узле сетки. Опциональное зеркалирование по плоскости YZ (флаг useSymmetry) позволяет хранить только половину карты для полей, антисимметричных по координате x.'),

  h2('2.2', 'Используемая информация'),
  tableCaption('2.1', 'Входные данные алгоритма'),
  table4col(
    ['Обозначение', 'Смысл', 'Ед.', 'Источник'],
    [
      ['p', 'точка опроса в локальной системе магнита', 'м', 'вызывающий код'],
      ['useSymmetry', 'флаг включения YZ-симметрии', '—', 'hardware_library.csv'],
      ['header.{minX..maxZ}', 'границы расчётной сетки', 'м', 'заголовок .bin'],
      ['header.{stepX..stepZ}', 'шаг сетки по осям', 'м', 'заголовок .bin'],
      ['header.{nx, ny, nz}', 'число узлов по осям', '—', 'заголовок .bin'],
      ['data[i·ny·nz + j·nz + k]', 'значение поля в узле (i,j,k)', 'Тл (В/м)', '.bin'],
    ]),

  h2('2.3', 'Результаты решения'),
  body('Возвращаемое значение — интерполированный вектор поля FieldVector в локальных координатах. Для точек, лежащих вне ограничивающего параллелепипеда сетки, возвращается нулевой вектор (0, 0, 0).'),

  h2('2.4', 'Математическое описание'),
  body('Поле задано в узлах регулярной декартовой сетки. Точка опроса p = (px, py, pz) в локальной системе попадает в единичную ячейку сетки; трилинейная интерполяция взвешивает значения в восьми углах этой ячейки.'),
  body('Дробные координаты ячейки. Сначала определяем ячейку, в которую попала точка:'),
  mathPara('i = ⌊', frac([sub('p', 'x'), ' − ', sub('x', 'min')], 'Δx'), '⌋,    ',
           sub('t', 'x'), ' = ',
           frac([sub('p', 'x'), ' − ', sub('x', 'min')], 'Δx'), ' − i'),
  mathPara('j = ⌊', frac([sub('p', 'y'), ' − ', sub('y', 'min')], 'Δy'), '⌋,    ',
           sub('t', 'y'), ' = ',
           frac([sub('p', 'y'), ' − ', sub('y', 'min')], 'Δy'), ' − j'),
  mathPara('k = ⌊', frac([sub('p', 'z'), ' − ', sub('z', 'min')], 'Δz'), '⌋,    ',
           sub('t', 'z'), ' = ',
           frac([sub('p', 'z'), ' − ', sub('z', 'min')], 'Δz'), ' − k'),
  body('Обозначения:'),
  dash('i, j, k — целые индексы нижнего угла ячейки (нумерация с нуля; верхний угол — (i+1, j+1, k+1)).'),
  dash('t_x, t_y, t_z ∈ [0, 1] — дробные смещения внутри ячейки; t_α = 0 — точка на нижней грани, t_α = 1 — на верхней.'),
  dash('x_min, y_min, z_min — координаты угла сетки с наименьшими индексами, м.'),
  dash('Δx, Δy, Δz — шаг сетки по каждой оси, м (равномерный по оси, но независимый между осями).'),
  body('Трилинейная интерполяция. Обозначим v_{abc} вектор поля в углу (i+a, j+b, k+c), где a, b, c ∈ {0, 1}. Тогда:'),
  mathPara('F(p) = ',
           new MathSum({
             children: [sub('w', 'abc'), '(', sub('t','x'),', ', sub('t','y'),', ', sub('t','z'),') · ', sub('v', 'abc')].map(_m),
             subScript: _arr('a,b,c ∈ {0,1}'),
           })),
  mathPara(sub('w', 'abc'), '(', sub('t','x'),', ', sub('t','y'),', ', sub('t','z'), ') = ',
           bra([sup(paren(['1 − ', sub('t','x')]), '1−a'), ' · ', sup(sub('t','x'), 'a')]),
           ' · ',
           bra([sup(paren(['1 − ', sub('t','y')]), '1−b'), ' · ', sup(sub('t','y'), 'b')]),
           ' · ',
           bra([sup(paren(['1 − ', sub('t','z')]), '1−c'), ' · ', sup(sub('t','z'), 'c')])),
  body('Обозначения:'),
  dash('v_{abc} — поле в одном из 8 углов (в Тл для B, в В/м для E).'),
  dash('w_{abc} — неотрицательный вес, линейно убывающий с расстоянием от точки опроса до соответствующего угла. Восемь весов всегда суммируются в единицу — результат является истинным взвешенным средним.'),
  body('В угловых точках (t_α ∈ {0, 1} по всем осям) ровно один вес равен 1, остальные — нулю, и F принимает значение в этом узле точно — интерполянт согласован в узлах.'),
  body('Правило YZ-симметрии. Если карта поля антисимметрична по x (хранится только половина объёма, x ∈ [0, x_max]), запрос с отрицательным x сворачивается в сохранённую половину:'),
  mathPara('F(', sub('p','x'),', ', sub('p','y'),', ', sub('p','z'), ') = ',
           sub('F', 'raw'), '(', sub('p','x'),', ', sub('p','y'),', ', sub('p','z'), '),    если ', sub('p','x'), ' ≥ 0'),
  mathPara('F(', sub('p','x'),', ', sub('p','y'),', ', sub('p','z'), ') = ',
           'diag(−1, +1, +1) · ', sub('F', 'raw'), '(−', sub('p','x'),', ', sub('p','y'),', ', sub('p','z'), '),    если ', sub('p','x'), ' < 0'),
  body('Обозначения:'),
  dash('F_raw — интерполяция, выполненная на сохранённой (px ≥ 0) половине сетки.'),
  dash('diag(−1, +1, +1) — диагональная матрица, инвертирующая x-компоненту и оставляющая остальные без изменений. Это отражает физическую антисимметрию B_x по плоскости YZ (тогда как B_y, B_z симметричны).'),

  h2('2.5', 'Алгоритм решения'),
  code(
    'Вход:  p, useSymmetry, grid (header + data)\n' +
    'Выход: вектор поля в точке p\n\n' +
    '1. mirrored ← false\n' +
    '2. если useSymmetry и p.x < 0:  p.x ← -p.x; mirrored ← true\n' +
    '3. если p вне [minX..maxZ]: вернуть {0,0,0}\n' +
    '4. fx ← (p.x - minX)/stepX;  i ← floor(fx);  tx ← fx - i\n' +
    '   fy ← (p.y - minY)/stepY;  j ← floor(fy);  ty ← fy - j\n' +
    '   fz ← (p.z - minZ)/stepZ;  k ← floor(fz);  tz ← fz - k\n' +
    '5. ограничить i, j, k интервалом [0, n-2] (граничный случай p == max)\n' +
    '6. извлечь 8 значений v000..v111 из data[]\n' +
    '7. покомпонентно вычислить F по формуле трилинейной интерполяции\n' +
    '8. если mirrored: F.x ← -F.x\n' +
    '9. вернуть F'),
];

// ---------------------------------------------------------------------------
// 3. Spatial bucket grid
// ---------------------------------------------------------------------------
const alg3 = [
  h1('3', 'Двумерный bucket-индекс магнитов'),
  body('Исходный код: src/LatticeManager.h:179 (buildSpatialIndex), src/LatticeManager.h:235 (getTotalFields).'),

  h2('3.1', 'Назначение и характеристика'),
  body('Ускорение поиска магнитов, потенциально вносящих вклад в поле в точке текущей позиции частицы. Без индекса симулятор должен был бы перебирать все ~2762 элемента кольца на каждом шаге. Индекс снижает число кандидатов до ~5 на шаг; на участках чистого дрейфа бакет пуст, и стоимость шага сводится к одному обращению к бакету без вычислений поля.'),
  body('Индекс строится один раз после завершения добавления элементов методом addElement(...) и далее используется в режиме «только чтение».'),

  h2('3.2', 'Используемая информация'),
  tableCaption('3.1', 'Входные данные алгоритма'),
  table3col(
    ['Обозначение', 'Смысл', 'Источник'],
    [
      ['elements[i].worldCenter', 'центр ограничивающей сферы магнита', 'computeBounds()'],
      ['elements[i].worldRadius', 'радиус ограничивающей сферы', 'computeBounds()'],
    ]),

  h2('3.3', 'Результаты решения'),
  tableCaption('3.2', 'Выходные данные алгоритма'),
  table2col(
    ['Обозначение', 'Смысл'],
    [
      ['bucketSize', 'длина ребра одного бакета (м)'],
      ['bxMin, byMin', 'нижний индекс сетки в единицах bucketSize'],
      ['bxN, byN', 'размеры сетки по осям x и y'],
      ['bucketGrid[bx·byN + by]', 'массив индексов элементов, пересекающих бакет'],
    ]),
  body('Запрос getTotalFields(pos) возвращает поле (B, E), просуммированное только по элементам того бакета, в котором находится точка pos.'),

  h2('3.4', 'Математическое описание'),
  body('У каждого магнита есть предвычисленная ограничивающая сфера с центром c_i (в мировых координатах) и радиусом r_i, полностью покрывающая его карту поля. Индекс работает в плоскости XY: кольцо лежит в z ≈ 0.'),
  body('Размер бакета. Пусть r_max — наибольший радиус ограничивающей сферы среди всех элементов. Шаг сетки:'),
  mathPara('h = max(1 м,  2 ', sub('r', 'max'), ')'),
  body('Обозначения:'),
  dash('h — длина ребра одного квадратного бакета, м.'),
  dash('r_max = max_i r_i — максимальный радиус по всем магнитам, м. Множитель 2 гарантирует, что XY-проекция любой ограничивающей окружности укладывается в группу 2×2 бакетов — магнит регистрируется не более чем в 4 бакетах.'),
  dash('Нижняя граница 1 м предотвращает патологически малые бакеты при наборе мелких магнитов (иначе массив бакетов разрастался бы).'),
  body('Координаты бакета. Мировой точке (x, y, z) сопоставляются индексы бакета (b_x, b_y):'),
  mathPara(sub('b', 'x'), ' = ⌊', frac('x', 'h'), '⌋ − ', subSup('b', 'x', 'min')),
  mathPara(sub('b', 'y'), ' = ⌊', frac('y', 'h'), '⌋ − ', subSup('b', 'y', 'min')),
  body('Обозначения:'),
  dash('floor(·) — округление вниз до ближайшего целого.'),
  dash('b_x^min, b_y^min — глобальные целочисленные смещения, подобранные при построении так, чтобы все занятые бакеты имели неотрицательные индексы. Вычисляются по самым малым x, y, достигаемым ограничивающими сферами.'),
  body('Магнит i регистрируется в каждом бакете (b_x, b_y), пересекающем диск { (x, y) : (x − c_{i,x})² + (y − c_{i,y})² ≤ r_i² }; по построению таких бакетов не более 2×2 = 4.'),

  h2('3.5', 'Алгоритм решения'),
  body('Построение индекса (выполняется однократно):'),
  code(
    '1. maxR ← max_i elements[i].worldRadius\n' +
    '2. bucketSize ← max(1.0, 2 · maxR)\n' +
    '3. xLo, xHi ← мин/макс по элементам от (worldCenter.x ± worldRadius)\n' +
    '   yLo, yHi ← аналогично по y\n' +
    '4. bxMin ← floor(xLo/bucketSize);  byMin ← floor(yLo/bucketSize)\n' +
    '   bxN ← floor(xHi/bucketSize) - bxMin + 1\n' +
    '   byN ← floor(yHi/bucketSize) - byMin + 1\n' +
    '5. bucketGrid ← массив пустых списков размера bxN·byN\n' +
    '6. для каждого элемента i:\n' +
    '       bxLo..bxHi ← диапазон бакетов по x\n' +
    '       byLo..byHi ← диапазон бакетов по y\n' +
    '       для bx ∈ [bxLo..bxHi], by ∈ [byLo..byHi]:\n' +
    '           bucketGrid[bx·byN + by].push_back(i)'),
  body('Запрос (на каждом шаге частицы):'),
  code(
    '1. bx ← floor(pos.x/bucketSize) - bxMin\n' +
    '   by ← floor(pos.y/bucketSize) - byMin\n' +
    '2. если (bx, by) вне сетки: B, E ← 0; вернуть\n' +
    '3. для каждого индекса i из bucketGrid[bx·byN + by]:\n' +
    '       elements[i].addContribution(pos, B, E)'),
  body('На участках дрейфа бакет пуст, поэтому внутренний цикл по кандидатам не выполняется и вычислений поля не происходит.'),
];

// ---------------------------------------------------------------------------
// 4. Magnet contribution
// ---------------------------------------------------------------------------
const alg4 = [
  h1('4', 'Вклад магнита с прескринингом по охватывающей сфере'),
  body('Исходный код: src/LatticeManager.h:114, метод MagneticElement::addContribution.'),

  h2('4.1', 'Назначение и характеристика'),
  body('Добавление вклада одного магнита в накапливающиеся суммы полей (B, E) в точке мирового опроса. Перед полным преобразованием координат и интерполяцией выполняется быстрая проверка по охватывающей сфере (7 арифметических операций), позволяющая отсечь заведомо далёкие элементы.'),

  h2('4.2', 'Используемая информация'),
  tableCaption('4.1', 'Входные данные алгоритма'),
  table3col(
    ['Обозначение', 'Смысл', 'Источник'],
    [
      ['globalPos', 'точка опроса в мировых координатах', 'цикл трекинга'],
      ['worldCenter, worldRadius²', 'параметры охватывающей сферы', 'computeBounds()'],
      ['transform', 'преобразование локальная↔мировая система', 'lattice_configs.csv'],
      ['map', 'карта поля магнита', 'hardware_library.csv'],
      ['useSymmetry', 'флаг YZ-симметрии (из колонки use_symmetry)', 'hardware_library.csv'],
      ['correctCoeff', 'калибровочный множитель', 'hardware_library + sp_elements'],
      ['fieldDirection', 'знак направления (±1)', 'hardware_library.csv'],
      ['currentScale', 'отношение текущего тока к номинальному', 'рассчитывается при загрузке'],
      ['isActive', 'признак активного магнита (применять currentScale)', 'hardware_library.csv'],
      ['type', 'MAGNETIC / ELECTRIC (из колонки field_type: B/M/MAGNETIC → MAGNETIC; E/ELECTRIC → ELECTRIC; отсутствие → MAGNETIC)', 'hardware_library.csv'],
    ]),
  body('addContribution маршрутизирует вклад по типу: MAGNETIC накапливается в totalB, ELECTRIC — в totalE. Все текущие устройства магнитные, поэтому totalE на практике остаётся нулём и шаги E-толчка Boris pusher\'а холостые, но логика поддерживает оба типа, и добавление электростатического элемента не потребует изменения кода.'),

  h2('4.3', 'Результаты решения'),
  body('Обновление по месту накопителя totalB (если type == MAGNETIC) или totalE (если type == ELECTRIC). Возвращаемого значения нет.'),

  h2('4.4', 'Математическое описание'),
  body('У каждого магнита собственная локальная система координат, в которой вычислена карта поля. Вклад магнита i в мировой точке x строится в четыре геометрических шага плюс масштабирование.'),
  body('Шаг 1 — прескрининг по охватывающей сфере. Магнит сразу отбрасывается, если точка опроса вне его сферы влияния:'),
  mathPara(sup(['|x − ', sub('c', 'i'), '|'], '2'), ' > ', sup(sub('r', 'i'), '2'),
           '    ⇒    вклад = 0'),
  dash('c_i — центр охватывающей сферы магнита i в мировых координатах, м.'),
  dash('r_i — радиус охватывающей сферы, м.'),
  body('Возведение в квадрат с обеих сторон позволяет обойтись без вычисления корня в горячем цикле.'),
  body('Шаг 2 — преобразование «мир → локал». Мировая точка переносится в локальную систему магнита:'),
  mathPara(sub('x', 'loc'), ' = ', sup(sub('R', 'i'), '⊤'), ' (x − ', sub('T', 'i'), ')'),
  dash('T_i — трансляция: мировая позиция локального начала координат магнита, м.'),
  dash('R_i — матрица поворота 3×3 (ортогональная, Rᵀ R = I), переводящая локальный вектор в мировой. Строится из углов yaw, pitch, roll из lattice_configs.csv.'),
  dash('Rᵀ — обратное преобразование (для ортогональной матрицы транспонирование совпадает с обратной), что отменяет размещение магнита.'),
  body('Шаг 3 — интерполяция поля в локальной системе (§ 2.4):'),
  mathPara(sub('F', 'loc'), ' = trilinear(', sub('map', 'i'), ', ', sub('x', 'loc'), ', ', sub('useSymmetry', 'i'), ')'),
  dash('map_i — карта поля магнита i (двоичная сетка векторов).'),
  dash('useSymmetry_i ∈ {true, false} — применять ли YZ-зеркалирование (см. § 2.4).'),
  body('Если F_loc = 0 (точка попала в нулевую зону сетки), возврат сразу.'),
  body('Шаг 4 — поворот вектора поля «локал → мир»:'),
  mathPara(sub('F', 'world'), ' = ', sub('R', 'i'), ' · ', sub('F', 'loc')),
  body('Трансляция здесь не применяется: магнитное поле — вектор, а не точка, и поворот — единственная действующая часть преобразования.'),
  body('Шаг 5 — масштабирование. Все индивидуальные калибровочные множители магнита перемножаются:'),
  mathPara(sub('scale', 'i'), ' = ',
           subSup('c', 'i', 'hw'), ' · ', subSup('c', 'i', 'sp'), ' · ', sub('d', 'i'),
           ' · '),
  mathPara('     · ',
           '{ ', frac(sub('I', 'i'), subSup('I', 'i', 'nom')), '   если магнит i активен;    1   иначе }'),
  dash('c_i^hw — устройственный коэффициент калибровки (hardware_library.correct_coeff), безразмерный. Одно значение на тип устройства.'),
  dash('c_i^sp — индивидуальный коэффициент калибровки (sp_elements.correct_coeff), безразмерный. По одному значению на пару (sp_type, sp_element_id); по умолчанию 1.0.'),
  dash('d_i ∈ {+1, −1} — флип направления (hardware_library.field_direction). Используется, когда физическая ориентация магнита отличается от заложенной в карту поля.'),
  dash('I_i — фактический ток через магнит, А (sp_elements.arg_val). Применяется только к активным электромагнитам; у постоянных магнитов тока нет и множитель равен 1.'),
  dash('I_i^nom — номинальный ток, при котором вычислялась карта поля, А (hardware_library.nominal_arg). Отношение I_i / I_i^nom масштабирует поле линейно по току (в ненасыщенном режиме).'),
  body('Шаг 6 — накопление в нужную сумму поля:'),
  mathPara(sub('B', 'total'), ' += ', sub('scale', 'i'), ' · ', sub('F', 'world'), ',    если ', sub('type', 'i'), ' = MAGNETIC'),
  mathPara(sub('E', 'total'), ' += ', sub('scale', 'i'), ' · ', sub('F', 'world'), ',    если ', sub('type', 'i'), ' = ELECTRIC'),
  dash('totalB, totalE — текущие суммы магнитного и электрического полей в точке опроса, Тл и В/м.'),
  dash('type_i — MAGNETIC (по умолчанию, все текущие устройства) или ELECTRIC (заготовка под будущие электростатические элементы), из hardware_library.field_type.'),

  h2('4.5', 'Алгоритм решения'),
  code(
    '1. dx,dy,dz ← globalPos − worldCenter\n' +
    '   если dx²+dy²+dz² > worldRadius²:  return\n' +
    '2. localPos ← transform.transformPointInv(globalPos)\n' +
    '3. localF   ← map.getInterpolatedField(localPos, useSymmetry)\n' +
    '4. если localF == 0:  return\n' +
    '5. globalF ← transform.rotateVector(localF)\n' +
    '6. scale ← correctCoeff · fieldDirection\n' +
    '   если isActive: scale *= currentScale\n' +
    '7. если type == MAGNETIC: totalB += scale · globalF\n' +
    '   иначе                  totalE += scale · globalF'),
];

// ---------------------------------------------------------------------------
// 5. Finish plane
// ---------------------------------------------------------------------------
const alg5 = [
  h1('5', 'Детектор пересечения финишной плоскости'),
  body('Исходный код: src/FinishPlane.h:51 (FinishPlane::hasCrossed), setup на строке 30.'),

  h2('5.1', 'Назначение и характеристика'),
  body('Детектирование момента, когда частица пересекает ориентированную плоскость, заданную тремя точками. Каждая частица обладает собственным состоянием детектора (счётчик вызовов, последнее знаковое расстояние), что обеспечивает корректную работу в многопоточном режиме OpenMP. Параметр skipSteps позволяет пропустить заданное число первых вызовов — это нужно, когда плоскость пересекает орбиту несколько раз (например, плоскость на половине кольца проходит и около старта, и около цели).'),

  h2('5.2', 'Используемая информация'),
  tableCaption('5.1', 'Входные данные алгоритма'),
  table3col(
    ['Обозначение', 'Смысл', 'Источник'],
    [
      ['a, b, c', 'три точки задания плоскости', 'experiments.csv'],
      ['startPos', 'начальная позиция частицы', 'particle_groups.csv'],
      ['skip', 'число пропускаемых вызовов', 'experiments.csv (max_turns)'],
      ['pos', 'текущая позиция частицы', 'цикл трекинга'],
    ]),

  h2('5.3', 'Результаты решения'),
  body('Метод hasCrossed(pos) возвращает true при первом вызове, в котором знаковое расстояние от точки pos до плоскости сменило знак относительно сохранённого lastD. Трекер использует этот флаг как сигнал останова интегрирования и записи траектории на диск.'),

  h2('5.4', 'Математическое описание'),
  body('Финишная плоскость задана тремя мировыми точками a, b, c, не лежащими на одной прямой. Они задают саму плоскость; алгоритм строит единичную нормаль и использует её для ориентированного «знакового» теста.'),
  body('Нормаль к плоскости. Из двух рёбер составляется векторное произведение:'),
  mathPara('n̂ = ', frac(['(b − a) × (c − a)'], ['|(b − a) × (c − a)|'])),
  body('Обозначения:'),
  dash('a, b, c — три задающие точки из experiments.csv, м.'),
  dash('n̂ — единичная нормаль к плоскости (модуль 1). Векторное произведение даёт вектор, ортогональный плоскости; знаменатель нормирует его длину.'),
  body('Соглашение о знаке. Ориентация n̂ из векторного произведения произвольна, поэтому она переворачивается так, чтобы частица начинала на «обратной» стороне:'),
  mathPara('n̂ ← −n̂,    если (', sub('x', '0'), ' − a) · n̂ > 0'),
  dash('x_0 — начальная позиция частицы, м (particle_groups.csv).'),
  dash('Скалярное произведение слева — знаковая проекция вектора x_0 − a на n̂. После этой коррекции частица по построению стартует со знакового расстояния меньше нуля.'),
  body('Знаковое расстояние произвольной точки x от плоскости:'),
  mathPara('d(x) = (x − a) · n̂'),
  dash('d > 0 — x находится на той же стороне, куда указывает n̂.'),
  dash('d < 0 — x на противоположной стороне.'),
  dash('|d| равно евклидову расстоянию от x до ближайшей точки плоскости.'),
  body('Детектирование пересечения. Между двумя последовательными отсчётами x_{n−1} и x_n (со знаковыми расстояниями d_{n−1} и d_n) пересечение фиксируется при смене знака:'),
  mathPara('( ', sub('d', 'n−1'), ' < 0   ∧   ', sub('d', 'n'), ' ≥ 0 )    ∨    ( ',
           sub('d', 'n−1'), ' ≥ 0   ∧   ', sub('d', 'n'), ' < 0 )'),
  body('Полуоткрытая форма (d_n ≥ 0 / d_n < 0) гарантирует, что вырожденный случай d_n = 0 не будет пропущен.'),

  h2('5.5', 'Алгоритм решения'),
  body('Инициализация (один раз на частицу):'),
  code(
    '1. u ← b − a;  v ← c − a;  n ← u × v\n' +
    '2. нормализовать n; при вырождении положить n ← (1, 0, 0)\n' +
    '3. если (startPos − a) · n > 0:  n ← −n\n' +
    '4. callCount ← 0;  lastD ← 0;  initialized ← false'),
  body('Запрос (на каждом шаге):'),
  code(
    '1. callCount ← callCount + 1\n' +
    '2. если callCount ≤ skipSteps:  вернуть false\n' +
    '3. d ← (pos − a) · n\n' +
    '4. если not initialized: lastD ← d; initialized ← true; вернуть false\n' +
    '5. crossed ← (lastD < 0 ∧ d ≥ 0) ∨ (lastD ≥ 0 ∧ d < 0)\n' +
    '6. lastD ← d\n' +
    '7. вернуть crossed'),
];

// ---------------------------------------------------------------------------
// 6. Distance to polyline
// ---------------------------------------------------------------------------
const alg6 = [
  h1('6', 'Расстояние от точки до полилинии проектной орбиты'),
  body('Исходный код: src/optimizer.cpp:66, функция distanceToPolyline.'),

  h2('6.1', 'Назначение и характеристика'),
  body('Вычисление минимального евклидова расстояния от точки до ломаной, заданной упорядоченным списком вершин. Используется в качестве амплитуды огибающей в envelope-компоненте функции потерь (см. раздел 7), а также в скрипте analyze_trajectory_stability.py для оценки отклонения траектории частицы от проектной орбиты.'),

  h2('6.2', 'Используемая информация'),
  tableCaption('6.1', 'Входные данные алгоритма'),
  table2col(
    ['Обозначение', 'Смысл'],
    [
      ['p', 'точка опроса в мировых координатах'],
      ['poly', 'упорядоченный список вершин полилинии (N ≥ 2)'],
    ]),

  h2('6.3', 'Результаты решения'),
  body('Скалярное расстояние ≥ 0. Если число вершин меньше 2, возвращается 0.'),

  h2('6.4', 'Математическое описание'),
  body('Полилиния — последовательность отрезков. Расстояние от точки до полилинии — минимум расстояний до отдельных отрезков; задача проекции на отрезок решается стандартно через клиппинг.'),
  body('Проекция на отрезок. Для отрезка [a, b] точка на бесконечной прямой, ближайшая к p, имеет параметр'),
  mathPara(sup('t', '*'), ' = ', frac('(p − a) · (b − a)', sup('|b − a|', '2'))),
  body('Обозначения:'),
  dash('a, b — концы отрезка, м.'),
  dash('p — точка опроса, м.'),
  dash('t* — параметр прямой; t* = 0 даёт a, t* = 1 даёт b, t* ∈ (0, 1) — точку между, t* вне [0, 1] — точку на продолжении.'),
  dash('|b − a|² — квадрат длины отрезка (пропускается, если ≤ 10⁻³⁰, чтобы корректно обработать вырожденный отрезок).'),
  body('Чтобы ограничиться самим отрезком (а не продолжением), t* клиппится:'),
  mathPara('t = clamp(', sup('t', '*'), ', 0, 1) = min(1,   max(0, ', sup('t', '*'), '))'),
  body('Ближайшая точка на отрезке к p:'),
  mathPara('c = a + t (b − a)'),
  mathPara(sub('d', 'i'), ' = |p − c|'),
  dash('c — основание перпендикулярной проекции (или ближайший конец, если перпендикуляр выходит за отрезок).'),
  dash('d_i — евклидово расстояние от p до отрезка i, м.'),
  body('Расстояние до полилинии — минимум по всем N − 1 отрезкам:'),
  mathPara('dist(p, Π) = ', sub('min', 'i = 0..N−2'), '   ', sub('d', 'i'),
           ' = ', rad([sub('min', 'i'), '  ', sup(['|p − ', sub('c', 'i'), '|'], '2')])),
  dash('Π = (p_0, p_1, …, p_{N−1}) — полилиния (упорядоченный список из N вершин).'),
  body('Корень берётся только один раз в конце; внутри цикла сравниваются квадраты расстояний — это и быстрее, и избавляет от ветвления для корня в горячем цикле.'),

  h2('6.5', 'Алгоритм решения'),
  code(
    'Вход:  p, poly[0..N-1]\n' +
    'Выход: расстояние dist ≥ 0\n\n' +
    '1. если N < 2:  вернуть 0\n' +
    '2. best² ← +∞\n' +
    '3. для i из [0..N-2]:\n' +
    '     a ← poly[i];  b ← poly[i+1]\n' +
    '     u ← b − a;  len² ← u·u\n' +
    '     t ← clamp( (p − a)·u / len²,  0,  1 ),  если len² > 1e-30,  иначе 0\n' +
    '     c ← a + t·u\n' +
    '     d² ← (p − c)·(p − c)\n' +
    '     если d² < best²: best² ← d²\n' +
    '4. вернуть √best²'),
];

// ---------------------------------------------------------------------------
// 7. Composite loss
// ---------------------------------------------------------------------------
const alg7 = [
  h1('7', 'Составная функция потерь оптимизатора'),
  body('Исходный код: src/optimizer.cpp, функция runTrackingAndLoss.'),

  h2('7.1', 'Назначение и характеристика'),
  body('Скалярная целевая функция, минимизируемая оптимизатором токов магнитов. Состоит из трёх аддитивных компонентов, отвечающих за различные физические дефекты пучка:'),
  dash('position_loss — частица должна проходить через центры магнитов;'),
  dash('direction_loss — пучок должен покидать решётку параллельно заданному направлению;'),
  dash('envelope_loss — бетатронные амплитуды не должны расти.'),
  body('Компонент position_loss обновляется инкрементально в цикле трекинга, без сохранения траектории. Direction- и envelope-компоненты вычисляются однократно на финальном шаге каждой частицы.'),

  h2('7.2', 'Используемая информация'),
  tableCaption('7.1', 'Входные данные алгоритма'),
  table3col(
    ['Обозначение', 'Смысл', 'Источник'],
    [
      ['траектории', '5 траекторий частиц (неявно через цикл)', 'трекинг'],
      ['magnet_sp_numbers', 'номера СП, участвующих в потере', 'JSON-конфиг'],
      ['exclude_device_types', 'исключаемые типы (ПМ-диполи, DQ и т. д.)', 'JSON-конфиг'],
      ['extra_target_magnets', 'виртуальные цели (x, y, z, weight)', 'JSON-конфиг'],
      ['direction_loss.reference', 'единичный эталонный вектор r̂', 'JSON-конфиг'],
      ['direction_loss.weight', 'вес w_dir', 'JSON-конфиг'],
      ['envelope_loss.weight', 'вес w_env', 'JSON-конфиг'],
      ['design_polyline', 'центры элементов + прелюдия/постлюдия', 'строится в начале запуска'],
    ]),

  h2('7.3', 'Результаты решения'),
  tableCaption('7.2', 'Выходные данные алгоритма'),
  table2col(
    ['Обозначение', 'Смысл'],
    [
      ['loss', 'скаляр для минимизации'],
      ['info.sum_d²', 'сумма квадратов расстояний для диагностики'],
      ['info.mean_mm, info.max_mm', 'среднее и максимальное отклонение в мм'],
      ['info.per_particle[i]', 'разбивка по частицам для history.json'],
    ]),

  h2('7.4', 'Математическое описание'),
  body('Полная целевая функция — сумма трёх слагаемых:'),
  mathPara('L = ', sub('L', 'pos'), ' + ', sub('L', 'dir'), ' + ', sub('L', 'env')),
  dash('L — скалярная потеря, минимизируемая оптимизатором (единицы — м² до взвешивания; веса нормируют слагаемые до сопоставимых порядков).'),
  body('1) Position loss — частицы должны проходить близко к центрам магнитов (и заданных пользователем виртуальных целей):'),
  mathPara(sub('L', 'pos'), ' = ',
           new MathSum({
             children: [
               new MathSum({
                 children: [sub('w', 'm'), '  ', sub('min', 't'), '  ',
                            sup(['|', sub('x', 'i'), '(t) − ', sub('c', 'm'), '|'], '2')].map(_m),
                 subScript: _arr('m ∈ M'),
               }),
             ],
             subScript: _arr('i = 1..N_p'),
           })),
  body('Обозначения:'),
  dash('N_p — число трекируемых частиц (5: опорная частица и 4 спутника).'),
  dash('M — множество «учитываемых магнитов» плюс виртуальных целей, после фильтрации по magnet_sp_numbers и exclude_device_types из JSON оптимизатора.'),
  dash('c_m — геометрический центр магнита m в мировых координатах, м.'),
  dash('w_m — индивидуальный вес цели (по умолчанию 1; targets из extra_target_magnets имеют свои бо́льшие веса, обычно 20 или 50, чтобы «пришпиливать» финальную сходимость).'),
  dash('x_i(t) — позиция частицы i на шаге интегрирования t, м.'),
  dash('min_t |x_i(t) − c_m|² — минимальный квадрат расстояния от траектории частицы i до магнита m по всей траектории. Считается инкрементально; хранится только текущий минимум, траектория не сохраняется.'),
  body('2) Direction loss — на финишной плоскости импульсы должны быть параллельны заданному эталонному направлению:'),
  mathPara(sub('L', 'dir'), ' = ', sub('w', 'dir'), '  ',
           new MathSum({
             children: [sup('sin', '2'), '(', sub('θ', 'i'), ')'].map(_m),
             subScript: _arr('i = 1..N_p'),
           })),
  mathPara('cos ', sub('θ', 'i'), ' = ',
           frac([subSup('p', 'i', 'final'), ' · r̂'], ['|', subSup('p', 'i', 'final'), '|'])),
  body('Обозначения:'),
  dash('w_dir — вес направленной компоненты (типично 10³); большой, потому что sin²(θ) безразмерно и обычно мало.'),
  dash('p_i_final — импульс частицы i в момент пересечения финишной плоскости, кг·м/с.'),
  dash('r̂ — единичный эталонный вектор из конфига оптимизатора (direction_loss.reference); обычно касательная к проектной орбите в точке финиша.'),
  dash('θ_i — угол между p_i_final и r̂, рад.'),
  dash('sin² θ_i = 1 − cos² θ_i — квадратично по θ при малых углах, что превращает компонент в least-squares угловую ошибку.'),
  body('3) Envelope loss — бетатронные осцилляции не должны расти при пролёте через решётку:'),
  mathPara(sub('L', 'env'), ' = ', sub('w', 'env'), '  ',
           new MathSum({
             children: [sup(['max(0,  ', subSup('A', 'i', 'final'), ' − ', subSup('A', 'i', 'init'), ')'], '2')].map(_m),
             subScript: _arr('i = 1..N_p'),
           })),
  mathPara(sub('A', 'i'), '(t) = dist( ', sub('x', 'i'), '(t), ', sub('Π', 'design'), ' )'),
  body('Обозначения:'),
  dash('w_env — вес огибающей, типично 10⁶ (большой, потому что A измеряется в метрах, а метровый рост недопустим).'),
  dash('A_i(t) — амплитуда огибающей частицы i — её расстояние от проектной полилинии в момент t, считается по алгоритму § 6.'),
  dash('Π_design — проектная полилиния (центры элементов решётки плюс прелюдийный/постлюдийный сегменты по касательной к орбите). Перестраивается в начале каждого прогона оптимизации.'),
  dash('A_i^init, A_i^final — амплитуды огибающей в начале и в конце прогона, м.'),
  dash('Обёртка max(0, ·) означает, что «затухание» (приближение пучка к орбите, A_i^final < A_i^init) НЕ штрафуется — штрафуется только рост.'),

  h2('7.5', 'Алгоритм решения'),
  body('А. Подготовка статических данных (выполняется один раз за вызов):'),
  code(
    '1. magnets_idx ← фильтрация элементов по (sp_number ∈ S) ∧\n' +
    '                                          (device_type ∉ exclude)\n' +
    '2. virtual_targets ← список (c, w) из extra_target_magnets\n' +
    '3. design_poly ← центры элементов 1..N + прелюдия\n' +
    '                  (экстраполяция первого сегмента на −2.94 м) +\n' +
    '                 постлюдия (экстраполяция последнего сегмента вперёд)'),
  body('Б. Цикл по каждой частице i:'),
  code(
    'minDist²[m] ← +∞ для каждого m из magnets_idx ∪ virtual_targets\n' +
    'A_init ← distanceToPolyline(x_init, design_poly)\n' +
    'для каждого шага до пересечения финишной плоскости:\n' +
    '    шаг Бориса (borisStep)\n' +
    '    для каждого m из бакета:\n' +
    '        d² ← |x − c_m|²\n' +
    '        если d² < minDist²[m]: minDist²[m] ← d²\n' +
    'A_final ← distanceToPolyline(x_final, design_poly)\n' +
    'θ_i     ← угол(p_final, r̂)'),
  body('В. Агрегация:'),
  code(
    'position_loss  = Σ_(i,m) w_m · minDist²[i, m]\n' +
    'direction_loss = w_dir · Σ_i sin²(θ_i)\n' +
    'envelope_loss  = w_env · Σ_i max(0, A_final − A_init)²\n' +
    'loss           = position_loss + direction_loss + envelope_loss'),
  body('Если хотя бы одна частица не достигает финишной плоскости за 50 миллионов шагов, функция возвращает признак сбоя, и оптимизатор отвергает текущее состояние.'),
];

// ---------------------------------------------------------------------------
// 8. Stochastic optimizer
// ---------------------------------------------------------------------------
const alg8 = [
  h1('8', 'Стохастический оптимизатор (hill-climber c киками)'),
  body('Исходный код: src/optimizer.cpp, функция main и вспомогательные функции.'),

  h2('8.1', 'Назначение и характеристика'),
  body('Однопоточный жадный hill-climber с периодическими кика́ми по всему вектору параметров для выхода из локальных минимумов. Оптимизатор изменяет токи магнитов (arg_val) и индивидуальные калибровочные множители (correct_coeff), перечисленные в parameter_groups_per_position JSON-конфига.'),
  body('Свойства реализации:'),
  dash('асимметричный шаг — на каждой итерации одна случайная группа возмущается на ±step_abs (ток) либо ±step_percent · current (множитель);'),
  dash('кики — каждые kick_every итераций все группы получают возмущение ±kick_scale одновременно;'),
  dash('жадный приём — новое состояние принимается, если loss уменьшился; иначе — откат;'),
  dash('контрольные точки — полная история сохраняется каждые CHECKPOINT_EVERY итераций для возможности восстановления после прерывания.'),

  h2('8.2', 'Используемая информация'),
  tableCaption('8.1', 'Входные данные алгоритма'),
  table3col(
    ['Обозначение', 'Смысл', 'Источник'],
    [
      ['experiment_id', 'идентификатор эксперимента', 'JSON-конфиг'],
      ['lattice_max_gid', 'трекинг только первых N gid решётки', 'JSON-конфиг'],
      ['n_iterations', 'общее число итераций', 'JSON-конфиг'],
      ['kick_every', 'период кика', 'JSON-конфиг'],
      ['kick_scale', 'амплитуда кика (А или доля)', 'JSON-конфиг'],
      ['parameter_groups_per_position', 'правила выбора параметров', 'JSON-конфиг'],
      ['step_abs / step_percent', 'шаг на группу', 'JSON-конфиг'],
      ['loss.*', 'веса и фильтры (см. п. 7.2)', 'JSON-конфиг'],
    ]),

  h2('8.3', 'Результаты решения'),
  tableCaption('8.2', 'Выходные данные алгоритма'),
  table2col(
    ['Файл', 'Содержание'],
    [
      ['sp_elements.csv', 'лучшее arg_val по каждой паре (sp_type, sp_eid)'],
      ['hardware_library.csv', 'лучшее correct_coeff по каждому device_type'],
      ['results/{name}.history.json', 'журнал итераций (loss, параметры)'],
      ['results/trajectories_exp{id}.bin', 'траектории при лучшем состоянии'],
    ]),

  h2('8.4', 'Математическое описание'),
  body('Оптимизатор обходит пространство параметров методом Монте-Карло с жадным приёмом.'),
  body('Обозначения:'),
  dash('θ ∈ ℝ^M — вектор всех изменяемых параметров (токи магнитов и индивидуальные коэффициенты), длина M.'),
  dash('G — число групп, объявленных в parameter_groups_per_position. Каждая группа g может объединять несколько слотов параметров; мутация группы g меняет все её слоты одновременно.'),
  dash('L(θ) — составная функция потерь из § 7.'),
  dash('θ_best — текущий лучший вектор параметров; L_best = L(θ_best).'),
  body('Обычная итерация. Случайно выбирается одна группа g; её параметры сдвигаются на индивидуальный шаг с случайным знаком:'),
  mathPara(sup(sub('θ', 'g'), '′'), ' = ', sub('θ', 'g'), ' + ', sub('ε', 'g'), ' · sign(u − ½)'),
  mathPara(sup(sub('θ', 'k'), '′'), ' = ', sub('θ', 'k'), ',    для всех k ≠ g'),
  body('Обозначения:'),
  dash('u ~ Uniform(0, 1) — свежее случайное число.'),
  dash('sign(u − ½) ∈ {−1, +1} — знак предложенного сдвига, плюс и минус равновероятны.'),
  dash('ε_g — размер шага для группы g, определяется типом колонки:'),
  mathPara(sub('ε', 'g'), ' = ', sub('step_abs', 'g'), ',    если группа g целит колонку arg_val'),
  mathPara(sub('ε', 'g'), ' = ', sub('step_percent', 'g'), '  · ', sub('θ', 'g'), ',    если группа g целит колонку correct_coeff'),
  dash('step_abs_g — абсолютный шаг в амперах (типично 0.1–1.0).'),
  dash('step_percent_g — относительный шаг (типично 0.05 ≡ 5%), пропорциональный текущему значению — относительное возмущение остаётся осмысленным при уменьшении θ_g.'),
  body('Кик-итерация. Каждые kick_every итераций ВСЕ группы возмущаются одновременно, чтобы выпрыгнуть из локальных минимумов:'),
  mathPara(sup(sub('θ', 'g'), '′'), ' = ', sub('θ', 'g'), ' + ', sub('δ', 'g'), ',    для всех g'),
  mathPara(sub('δ', 'g'), ' ~ Uniform(−Δ, +Δ)'),
  body('Обозначения:'),
  dash('Δ = kick_scale — амплитуда кика (А или доля, в зависимости от типа колонки группы).'),
  dash('Независимые равномерные выборки означают, что каждая группа движется независимо; совместное предложение достаточно крупное, чтобы покинуть бассейн, который пошаговый поиск никогда не покинет.'),
  body('Жадный приём. Вне зависимости от типа итерации, предложение принимается только при уменьшении лучшего значения потерь:'),
  mathPara('θ,  ', sub('L', 'best'), '  ←  θ′,  L(θ′),    если L(θ′) < ', sub('L', 'best')),
  mathPara('θ,  ', sub('L', 'best'), '  ←  θ,  ', sub('L', 'best'), ',    иначе (откат)'),
  body('Это чисто hill-climber — без «температуры», без Metropolis-приёма худших состояний. Кики и дают шанс покинуть локальные минимумы.'),

  h2('8.5', 'Алгоритм решения'),
  code(
    '1. Загрузить конфигурацию:\n' +
    '     - прочитать JSON, заполнить OptimConfig\n' +
    '     - прочитать sp_elements, hardware_library, lattice_configs,\n' +
    '       particle_groups\n' +
    '     - отфильтровать строки решётки по lattice_id и lattice_max_gid\n' +
    '2. Построить ускоритель:\n' +
    '     - загрузить карты поля в общий кеш (один .bin на device_type)\n' +
    '     - для каждой строки: создать MagneticElement, computeBounds()\n' +
    '     - acc.buildSpatialIndex()\n' +
    '3. Построить список параметров:\n' +
    '     - для каждой записи parameter_groups_per_position:\n' +
    '         найти соответствующие (sp_type, sp_eid)\n' +
    '         собрать индексы элементов\n' +
    '         сохранить ссылки на изменяемые поля\n' +
    '4. Начальная оценка:\n' +
    '     L_best ← runTrackingAndLoss(θ)\n' +
    '     при сбое — прервать\n' +
    '5. Цикл i = 1..n_iterations:\n' +
    '     если i % kick_every == 0:\n' +
    '         для каждой группы g: θ_g += Uniform(−kick_scale, +kick_scale)\n' +
    '     иначе:\n' +
    '         выбрать случайную группу g\n' +
    '         θ_g_old ← θ_g\n' +
    '         θ_g ← θ_g + step · sign(rand − 0.5)\n' +
    '     L_new ← runTrackingAndLoss(θ)\n' +
    '     если L_new — сбой или L_new ≥ L_best:\n' +
    '         откатить мутации этой итерации\n' +
    '     иначе:\n' +
    '         L_best ← L_new\n' +
    '         добавить (i, L_best, snapshot θ) в историю\n' +
    '     если i % CHECKPOINT_EVERY == 0:\n' +
    '         сохранить history.json (атомарно через временный файл)\n' +
    '6. Финальная фаза:\n' +
    '     - записать лучшее θ обратно в sp_elements.csv и\n' +
    '       hardware_library.csv\n' +
    '     - пересобрать lattice_configs.csv (сохранить суб-решётки)\n' +
    '     - выполнить финальный прогон, записать trajectories_exp{id}.bin\n' +
    '     - сохранить итоговый history.json'),
];

// ---------------------------------------------------------------------------
// 9. LTE -> lattice
// ---------------------------------------------------------------------------
const alg9 = [
  h1('9', 'Генерация решётки из LTE-описания'),
  body('Исходный код: lte_to_lattice.py.'),

  h2('9.1', 'Назначение и характеристика'),
  body('Разбор файлов формата Elegant LTE (SYLA_desc.lte, SYLA_line.lte), классификация элементов согласно конфигурации lte_mapping.json, развёртка кольца вдоль хорды и формирование редактируемых шаблонов sp_elements.csv и геометрической развёртки lattice_layout.csv.'),
  body('Особенности:'),
  dash('секционирование DL/JL — длинные изгибающие магниты делятся на 5 фрагментов, размещаемых на хорде вместо дуги;'),
  dash('сохранение пользовательских arg_val — при повторных запусках значения, ранее введённые в sp_elements.csv, сохраняются.'),

  h2('9.2', 'Используемая информация'),
  tableCaption('9.1', 'Входные данные алгоритма'),
  table2col(
    ['Путь', 'Содержание'],
    [
      ['workspace/lte/SYLA_desc.lte', 'описания устройств (тип, длина, K1, углы)'],
      ['workspace/lte/SYLA_line.lte', 'последовательность элементов кольца'],
      ['workspace/configs/lte_mapping.json', 'соответствие LTE-типа → device_type_id, регулярные выражения, пороги классификации'],
    ]),

  h2('9.3', 'Результаты решения'),
  tableCaption('9.2', 'Выходные данные алгоритма'),
  table2col(
    ['Путь', 'Содержание'],
    [
      ['workspace/tables/sp_elements.csv', 'шаблоны (sp_type, sp_element_id, device_type_id, arg_val, correct_coeff)'],
      ['workspace/tables/lattice_layout.csv', 'развёртка (global_id, sp_number, sp_type, sp_element_id, x, y, z, yaw, pitch, roll)'],
    ]),

  h2('9.4', 'Математическое описание'),
  body('Длина хорды поворотного магнита. У поворотного магнита заданы длина дуги L и угол поворота α, радиус кривизны ρ = L/α. Прямая хорда, соединяющая входную и выходную грани дуги, имеет длину'),
  mathPara(sub('L', 'chord'), ' = 2 ρ sin(α/2) = ',
           frac('2 L sin(α/2)', 'α')),
  body('Обозначения:'),
  dash('L — длина дуги поворотного магнита, м (из параметра L= LTE).'),
  dash('α — полный угол поворота, рад (из параметра ANGLE= LTE). α > 0 — поворот влево, α < 0 — вправо.'),
  dash('ρ = L/α — радиус кривизны, м.'),
  dash('L_chord — длина хорды, по которой выпрямляется магнит (чуть короче L при α ≠ 0). При α → 0 хорда стремится к L, как и ожидается.'),
  body('Магнит заменяется 5 тонкими ломтиками, равномерно расположенными по хорде в параметрах τ ∈ {0, 0.25, 0.5, 0.75, 1.0}, каждый ориентирован по направлению хорды.'),
  body('Обход кольца. Развёртка идёт сегмент за сегментом. Каждый LTE-элемент длины ℓ даёт в плоскости поворот Δψ (0 для дрейфа, α для поворотного магнита) и сдвиг геометрической середины. Обновление состояния:'),
  mathPara('ψ′ = ψ + Δψ'),
  mathPara('x′ = x + ℓ cos( ψ + ', frac('Δψ', '2'), ' )'),
  mathPara('y′ = y + ℓ sin( ψ + ', frac('Δψ', '2'), ' )'),
  body('Обозначения:'),
  dash('ψ — yaw-угол (поворот вокруг мировой оси z), рад.'),
  dash('Δψ — yaw-поворот, наработанный этим элементом, рад (равен α для поворотного магнита, 0 для прямого элемента).'),
  dash('ℓ — длина элемента по хорде (для поворотного магнита это пятая часть L_chord на ломтик).'),
  dash('x, y — накопленные мировые координаты текущей позиции, м.'),
  body('Использование «среднего» yaw ψ + Δψ/2 в тригонометрических аргументах (вместо входного ψ или выходного ψ + Δψ) сохраняет полилинию касательной к истинной дуге на концах и симметрично уравновешивает ошибку кривизны вдоль ломтика.'),

  h2('9.5', 'Алгоритм решения'),
  code(
    '1. Чтение JSON-конфига → mapping, regex-шаблоны, пороги\n' +
    '2. Разбор SYLA_desc.lte:\n' +
    '     для каждой строки "ИМЯ: ТИП, key1=val1, ..."\n' +
    '         классифицировать ТИП через конфиг (DL/JL/SBM/DP/DQ/QP/SQ/COR/MARK/...)\n' +
    '         сохранить (ИМЯ, классиф_тип, длина, угол, K1, ...)\n' +
    '3. Разбор SYLA_line.lte:\n' +
    '     развернуть вложенные LINE-выражения в плоскую\n' +
    '     последовательность (ИМЯ, повтор)\n' +
    '4. Классификация в (sp_type, sp_element_id):\n' +
    '     standard — обычные СП, повторяющиеся 1..N\n' +
    '     junction — стыковочные СП 18-21 (отличающаяся геометрия)\n' +
    '5. Развёртка кольца:\n' +
    '     для каждого элемента в плоской последовательности:\n' +
    '         если DL/JL: выдать 5 фрагментов на хорде\n' +
    '         иначе:      выдать 1 фрагмент\n' +
    '         обновить (x, y, yaw) по формулам п. 9.4\n' +
    '6. Прочитать существующий sp_elements.csv (если есть):\n' +
    '     сохранить arg_val для совпадающих ключей (sp_type, sp_eid)\n' +
    '7. Записать sp_elements.csv (шаблоны, по одной строке на ключ)\n' +
    '8. Записать lattice_layout.csv (полная развёртка)'),
];

// ---------------------------------------------------------------------------
// 10. Sub-lattice assembly
// ---------------------------------------------------------------------------
const alg10 = [
  h1('10', 'Сборка таблицы конфигурации и суб-решёток'),
  body('Исходный код: assemble_lattice_configs.py.'),

  h2('10.1', 'Назначение и характеристика'),
  body('Объединение редактируемого шаблона sp_elements.csv с геометрической развёрткой lattice_layout.csv в итоговую входную таблицу симулятора lattice_configs.csv. Дополнительно генерирует укороченные копии решётки (id 102, 103, 104, ...) на основании sub_lattices.json, содержащие только первые max_gid строк основной решётки. Это позволяет запускать эксперименты на части кольца без повторной развёртки LTE.'),

  h2('10.2', 'Используемая информация'),
  tableCaption('10.1', 'Входные данные алгоритма'),
  table2col(
    ['Путь', 'Содержание'],
    [
      ['workspace/tables/sp_elements.csv', 'шаблоны СП (редактируется вручную)'],
      ['workspace/tables/lattice_layout.csv', 'геометрия из lte_to_lattice.py'],
      ['workspace/configs/lte_mapping.json', 'lattice_id основной решётки'],
      ['workspace/configs/sub_lattices.json', 'список {lattice_id, max_gid, description}'],
    ]),

  h2('10.3', 'Результаты решения'),
  tableCaption('10.2', 'Выходные данные алгоритма'),
  table2col(
    ['Путь', 'Содержание'],
    [
      ['workspace/tables/lattice_configs.csv',
       'по одной строке на пару (lattice, instance), колонки: lattice_id, instance_id, device_type_id, x, y, z, yaw, pitch, roll, arg_val, correct_coeff'],
    ]),

  h2('10.4', 'Математическое описание'),
  body('Арифметика отсутствует — выполняется реляционное соединение (join) lattice_layout (геометрия) и sp_elements (токи и калибровки), плюс объединение с укороченными копиями основной решётки.'),
  body('В нотации построения множества:'),
  code(
    'Output = R_master  ∪  ⋃_{s ∈ SubLattices} R_s\n\n' +
    'R_master = { ( L_0, ℓ.global_id, s.device_type_id, ℓ.q,\n' +
    '              s.arg_val, s.correct_coeff )\n' +
    '             | ℓ ∈ lattice_layout,\n' +
    '               s ∈ sp_elements,\n' +
    '               (s.sp_type, s.sp_element_id) =\n' +
    '                  (ℓ.sp_type, ℓ.sp_element_id) }\n\n' +
    'R_s = { ( s.lattice_id, остальные поля строки )\n' +
    '        | r ∈ R_master,  r.instance_id ≤ s.max_gid }'),
  body('Обозначения:'),
  dash('L_0 = 101 — идентификатор основной решётки (из lte_mapping.json).'),
  dash('ℓ — строка lattice_layout.csv.'),
  dash('ℓ.q = (ℓ.x, ℓ.y, ℓ.z, ℓ.yaw, ℓ.pitch, ℓ.roll) — геометрический кортеж элемента.'),
  dash('s — в первом блоке: строка sp_elements.csv; во втором: дескриптор суб-решётки из sub_lattices.json.'),
  dash('Условие соединения — по составному ключу (sp_type, sp_element_id): каждая строка layout берёт свой device_type_id, ток и коэффициент из соответствующего SP-шаблона.'),
  dash('Объединение по суб-решёткам создаёт дополнительные строки с новым lattice_id, но только для instance_id ≤ max_gid, то есть с префиксом основного кольца.'),

  h2('10.5', 'Алгоритм решения'),
  code(
    '1. sp_dict ← load_sp_elements()              ; ключ (sp_type, sp_eid)\n' +
    '2. output_rows ← []\n' +
    '3. для каждой строки r из lattice_layout.csv:\n' +
    '     key ← (r.sp_type, r.sp_element_id)\n' +
    '     если key не в sp_dict: выдать warning и пропустить\n' +
    '     dev, arg, coeff ← sp_dict[key]\n' +
    '     добавить запись {LATTICE_ID, r.global_id, dev, r.x..roll, arg, coeff}\n' +
    '4. для каждого (sub_id, max_gid) из SUB_LATTICES:\n' +
    '     для r из output_rows, где r.instance_id ≤ max_gid:\n' +
    '         добавить копию с lattice_id ← sub_id\n' +
    '5. Запись CSV с повторными попытками при PermissionError\n' +
    '   (блокировка файла в Windows)'),
];

// ---------------------------------------------------------------------------
// 11. Envelope + FFT
// ---------------------------------------------------------------------------
const alg11 = [
  h1('11', 'Анализ огибающей и бетатронных частот по FFT'),
  body('Исходный код: analyze_trajectory_stability.py.'),

  h2('11.1', 'Назначение и характеристика'),
  body('Постобработка бинарного файла траектории для оценки динамической устойчивости пучка:'),
  dash('построение огибающей A_i(t) = dist(x_i(t), design_polyline) для каждой из 5 частиц;'),
  dash('подгонка линейной и экспоненциальной моделей роста A_i(t) с оценкой коэффициента детерминации R²;'),
  dash('FFT-анализ A_i(t) для определения частоты бетатронных колебаний (tune);'),
  dash('выработка вердикта по каждой частице: STABLE / OSCILLATING / GROWING / DRIFTED.'),

  h2('11.2', 'Используемая информация'),
  tableCaption('11.1', 'Входные данные алгоритма'),
  table2col(
    ['Путь', 'Содержание'],
    [
      ['workspace/results/trajectories_exp{id}.bin', '5 траекторий частиц'],
      ['workspace/tables/lattice_configs.csv', 'центры элементов → проектная полилиния'],
      ['workspace/tables/experiments.csv', 'идентификатор использованной решётки'],
    ]),

  h2('11.3', 'Результаты решения'),
  body('Для каждой частицы формируются:'),
  dash('график огибающей A(t) в matplotlib;'),
  dash('коэффициенты линейной подгонки: наклон m_lin и свободный член A₀_lin;'),
  dash('коэффициенты экспоненциальной подгонки: показатель γ_exp и амплитуда A₀_exp вместе с оценкой R²;'),
  dash('FFT-спектр и доминирующая частота f_betatron;'),
  dash('классификационный вердикт.'),

  h2('11.4', 'Математическое описание'),
  body('Амплитуда огибающей. Для каждой частицы i и каждого сохранённого отсчёта траектории k:'),
  mathPara(sub('A', 'i'), '(', sub('t', 'k'), ') = dist( ', sub('x', 'i'), '(', sub('t', 'k'), '),  ', sub('Π', 'design'), ' )'),
  body('Обозначения:'),
  dash('x_i(t_k) — позиция частицы i в k-й сохранённой временной точке, м.'),
  dash('Π_design — проектная полилиния (центры элементов + прелюдия/постлюдия по касательной).'),
  dash('A_i(t_k) — мгновенное евклидово расстояние от траектории до проектной орбиты, м. Считается по алгоритму § 6.'),
  body('Линейная подгонка (базовая дрейф-модель):'),
  mathPara('A(t) ≈ m t + ', sub('A', '0')),
  mathPara('(m, ', sub('A', '0'), ') = argmin   ',
           new MathSum({
             children: [sup(paren(['A(', sub('t', 'k'), ') − m  ', sub('t', 'k'), ' − ', sub('A', '0')]), '2')].map(_m),
             subScript: _arr('k'),
           })),
  body('Обозначения:'),
  dash('m — линейный темп роста огибающей, м/с. Положительный m означает дрейф пучка от орбиты.'),
  dash('A_0 — экстраполированная огибающая в t = 0, м.'),
  dash('Минимизация — стандартный метод наименьших квадратов (замкнутая формула через numpy.polyfit(..., deg=1)).'),
  body('Экспоненциальная подгонка (модель роста/затухания):'),
  mathPara('A(t) ≈ ', sub('A', '0'), ' · ', sup('e', 'γ t')),
  mathPara('ln A(t) ≈ ln ', sub('A', '0'), ' + γ t'),
  body('Обозначения:'),
  dash('γ — экспоненциальный темп роста, обратные секунды (своего рода Ляпунов-подобный показатель огибающей). γ > 0 — амплитуда растёт; γ < 0 — затухает.'),
  dash('Подгонка линейна по ln A, поэтому выполняется в лог-пространстве методом наименьших квадратов по точкам с A > A_floor (чтобы избежать log 0).'),
  body('Качество экспоненциальной подгонки оценивается коэффициентом детерминации:'),
  mathPara(sup('R', '2'), ' = 1 − ', frac(sub('SS', 'res'), sub('SS', 'tot'))),
  mathPara(sub('SS', 'res'), ' = ',
           new MathSum({
             children: [sup(paren(['ln A(', sub('t', 'k'), ') − ', sub('lnA', 'fit'), '(', sub('t', 'k'), ')']), '2')].map(_m),
             subScript: _arr('k'),
           })),
  mathPara(sub('SS', 'tot'), ' = ',
           new MathSum({
             children: [sup(paren(['ln A(', sub('t', 'k'), ') − ⟨ln A⟩']), '2')].map(_m),
             subScript: _arr('k'),
           })),
  body('Обозначения:'),
  dash('SS_res — сумма квадратов остатков (ошибка подгонки).'),
  dash('SS_tot — суммарная сумма квадратов (вариация данных).'),
  dash('⟨ln A⟩ — среднее значение логарифма амплитуды.'),
  dash('R² ∈ (−∞, 1] — ближе к 1 = лучшая подгонка; R² < 0 — подгонка хуже, чем использование среднего.'),
  body('Низкий R² дисквалифицирует оценку γ, даже если она численно большая.'),
  body('Бетатронная частота через FFT. Дискретное преобразование Фурье огибающей при равномерном шаге Δt_sample = t_k − t_{k−1}:'),
  mathPara('Â(', sub('f', 'n'), ') = ',
           new MathSum({
             children: ['A(', sub('t', 'k'), ') · ', sup('e', ['−2π i  ', sub('f', 'n'), '  ', sub('t', 'k')])].map(_m),
             subScript: _arr('k = 0..K−1'),
           })),
  mathPara(sub('f', 'n'), ' = ', frac('n', ['K  ', sub('Δt', 'sample')]), ',    n = 0, 1, …, K−1'),
  body('Обозначения:'),
  dash('K — число отсчётов огибающей в траектории.'),
  dash('f_n — дискретный частотный бин, Гц.'),
  dash('Â(f_n) — комплексный коэффициент Фурье на частоте f_n.'),
  body('Доминирующая бетатронная частота — бин с максимальной амплитудой, исключая DC-бин (он несёт только среднее):'),
  mathPara(sub('f', 'betatron'), ' = ', sub('argmax', [sub('f', 'n'), ' > ', sub('f', 'min')]),
           '   |Â(', sub('f', 'n'), ')|'),
  body('Обозначения:'),
  dash('f_min — малая пороговая частота, маскирующая DC-компоненту (и, опционально, несколько соседних бинов, чтобы избавиться от утечки от ненулевого среднего).'),
  dash('f_betatron — частота сильнейшего спектрального пика, Гц. Для перевода в дробный «тюн» (tune) нужно поделить на частоту обращения кольца; анализатор делает это, когда длина решётки известна.'),

  h2('11.5', 'Алгоритм решения'),
  code(
    '1. Загрузить траектории (5 массивов точек (Ni, 3)).\n' +
    '2. Определить lattice_id по experiments.csv, загрузить lattice_configs,\n' +
    '   построить design_poly (центры + прелюдия/постлюдия).\n' +
    '3. Для каждой частицы i:\n' +
    '     a. A_i ← [distanceToPolyline(p, design_poly) для p из traj_i]\n' +
    '     b. t_i ← равномерный шаг (dt × traj_stride × stride)\n' +
    '     c. линейная подгонка: m, A₀ ← polyfit(t_i, A_i, deg=1)\n' +
    '     d. экспоненциальная подгонка (только точки с A > A_floor):\n' +
    '            log A = γ·t + log A₀; вычислить R²\n' +
    '     e. FFT(A_i); маскировать DC-бин; выбрать максимум как f_betatron\n' +
    '     f. классификация:\n' +
    '            DRIFTED      если A_init ≈ 0 и m_lin > δ\n' +
    '            GROWING      если γ > γ_threshold и R² ≥ R²_min\n' +
    '            OSCILLATING  если есть выраженный пик FFT и γ ≤ γ_threshold\n' +
    '            STABLE       иначе\n' +
    '4. Построить графики, вывести таблицу вердиктов.'),
];

// ---------------------------------------------------------------------------
// Приложения
// ---------------------------------------------------------------------------
const appA = [
  h1Unnumbered('Приложение А. Файлы конфигурации'),
  body('В таблице А.1 приведены конфигурационные файлы проекта и их назначение.'),
  tableCaption('А.1', 'Файлы конфигурации'),
  table3col(
    ['Файл', 'Используется', 'Назначение'],
    [
      ['workspace/configs/lte_mapping.json', 'lte_to_lattice.py, assemble_lattice_configs.py',
        'Соответствие LTE-типов и lattice_id'],
      ['workspace/configs/sub_lattices.json', 'assemble_lattice_configs.py, src/optimizer.cpp',
        'Определения суб-решёток'],
      ['workspace/configs/visualizer.json', 'visualize.py, live_visualizer.py',
        'Камера, цвета, прозрачность, проекция'],
      ['workspace/optim_inputs/*.json', 'prepare_optim_config.py → SynchrotronOptimizer[.exe]',
        'Настройки оптимизатора для эксперимента'],
      ['_common.py', 'все Python-скрипты верхнего уровня',
        'Кросс-платформенные хелперы: суффикс .exe, корень проекта, выбор matplotlib-бэкенда'],
    ]),
];

const appB = [
  h1Unnumbered('Приложение Б. Файлы данных'),
  body('В таблице Б.1 приведены файлы данных и их форматы.'),
  tableCaption('Б.1', 'Файлы данных'),
  table3col(
    ['Файл', 'Формат', 'Описание'],
    [
      ['workspace/fields_bin/*.bin', 'двоичный: GridHeader + FieldVector[nx·ny·nz]',
        'карта поля одного магнита'],
      ['workspace/tables/hardware_library.csv', 'CSV «;», 10 колонок',
        'калибровка по device_type'],
      ['workspace/tables/sp_elements.csv', 'CSV «;»',
        'шаблон СП (редактируется)'],
      ['workspace/tables/lattice_layout.csv', 'CSV «;»',
        'развёртка кольца (генерируется)'],
      ['workspace/tables/lattice_configs.csv', 'CSV «;»',
        'входная таблица симулятора (генерируется)'],
      ['workspace/tables/experiments.csv', 'CSV «;»',
        'определения экспериментов'],
      ['workspace/tables/particle_groups.csv', 'CSV «;»',
        'начальные условия частиц'],
      ['workspace/results/trajectories_exp{id}.bin', 'двоичный',
        'выход SynchrotronTracker[.exe]'],
      ['workspace/results/{name}.history.json', 'JSON',
        'выход SynchrotronOptimizer[.exe]'],
    ]),
  body('Колонки hardware_library.csv (в порядке файла): device_type_id, Active, field_type, map_raw_path, map_bin_path, stl_path, use_symmetry, nominal_arg, field_direction, correct_coeff. Флаги field_type (B/E) и use_symmetry (1/0) читаются в MagneticElement::type и MagneticElement::useSymmetry при загрузке. Все текущие строки имеют B / 1, но маршрутизация в addContribution и зеркальная ветка в FieldMap::getInterpolatedField реально учитывают оба флага — для добавления электростатического элемента или асимметричной карты поля изменения кода не требуется.'),
];

// ---------------------------------------------------------------------------
// Содержание (статическое)
// ---------------------------------------------------------------------------
const toc = [
  h1Unnumbered('Содержание', { pageBreakBefore: false }),
  ...TOC_ENTRIES.map(([num, text]) => tocItem(num, text, '')),
];

// ---------------------------------------------------------------------------
// Сборка документа
// ---------------------------------------------------------------------------

// Колонтитулы: на титуле — пустой; на остальных — номер по центру внизу.
const footerWithPageNumber = new Footer({
  children: [ new Paragraph({
    alignment: AlignmentType.CENTER,
    indent: { firstLine: 0 },
    children: [ new TextRun({ children: [PageNumber.CURRENT],
                              font: TIMES, size: SIZE_BODY }) ],
  }) ],
});

const titleSection = {
  properties: {
    page: {
      size:   { width: PAGE.width, height: PAGE.height },
      margin: { top: PAGE.marginTop, right: PAGE.marginRight,
                bottom: PAGE.marginBottom, left: PAGE.marginLeft },
    },
    titlePage: true,
  },
  headers: { default: new Header({ children: [new Paragraph('')] }) },
  footers: {
    first:   new Footer({ children: [new Paragraph('')] }),  // титул без номера
    default: footerWithPageNumber,
  },
  children: titlePage,
};

const mainSection = {
  properties: {
    page: {
      size:   { width: PAGE.width, height: PAGE.height },
      margin: { top: PAGE.marginTop, right: PAGE.marginRight,
                bottom: PAGE.marginBottom, left: PAGE.marginLeft },
      pageNumbers: { start: 2 },
    },
  },
  footers: { default: footerWithPageNumber },
  children: [
    ...toc,
    ...intro,
    ...alg1,  ...alg2,  ...alg3, ...alg4,  ...alg5,
    ...alg6,  ...alg7,  ...alg8, ...alg9,  ...alg10, ...alg11,
    ...appA, ...appB,
  ],
};

const doc = new Document({
  creator: 'SynchrotronSim',
  title:   'Описание алгоритмов',
  description: 'Описание алгоритмов SynchrotronTracker / SynchrotronOptimizer по РД 50-34.698-90',
  styles: {
    default: { document: { run: { font: TIMES, size: SIZE_BODY } } },
    paragraphStyles: [
      { id: 'Heading1', name: 'Heading 1', basedOn: 'Normal', next: 'Normal',
        quickFormat: true,
        run: { size: SIZE_BODY, bold: true, font: TIMES },
        paragraph: { spacing: { before: 480, after: 240, line: 360, lineRule: 'auto' },
                     outlineLevel: 0, alignment: AlignmentType.CENTER } },
      { id: 'Heading2', name: 'Heading 2', basedOn: 'Normal', next: 'Normal',
        quickFormat: true,
        run: { size: SIZE_BODY, bold: true, font: TIMES },
        paragraph: { spacing: { before: 280, after: 140, line: 360, lineRule: 'auto' },
                     outlineLevel: 1 } },
    ],
  },
  sections: [titleSection, mainSection],
});

const outPath = path.join(__dirname, 'algorithms.docx');
Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(outPath, buf);
  console.log(`Wrote ${outPath} (${(buf.length/1024).toFixed(1)} KiB)`);
}).catch(e => {
  console.error('Error generating docx:', e);
  process.exit(1);
});
