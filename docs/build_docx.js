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
  Bookmark, InternalHyperlink,
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

// Авто-нумерация формул в формате (X.Y), где X — номер главы (задаётся
// через startSection), Y — порядковый номер формулы внутри главы.
let _currentSection = null;
const _formulaCounters = new Map();
function startSection(num) {
  _currentSection = num;
  if (num !== null && !_formulaCounters.has(num)) _formulaCounters.set(num, 0);
}
function _nextFormulaLabel() {
  if (_currentSection === null) return null;
  const n = (_formulaCounters.get(_currentSection) || 0) + 1;
  _formulaCounters.set(_currentSection, n);
  return `(${_currentSection}.${n})`;
}

// Display-формула. Если активна нумерация (startSection вызван), формула
// получает метку (X.Y), прижатую к правому краю содержимой области.
function mathPara(...nodes) {
  const label = _nextFormulaLabel();
  const mathNode = new MathBlock({ children: nodes.map(_m) });
  if (label === null) {
    // Без номера — центрировано.
    return new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { line: 360, lineRule: 'auto', before: 160, after: 160 },
      indent: { firstLine: 0 },
      children: [ mathNode ],
    });
  }
  // С номером: левое выравнивание + центрирующий tab + правый tab для (X.Y).
  return new Paragraph({
    alignment: AlignmentType.LEFT,
    spacing: { line: 360, lineRule: 'auto', before: 160, after: 160 },
    indent: { firstLine: 0 },
    tabStops: [
      { type: TabStopType.CENTER, position: Math.round(CONTENT_W / 2) },
      { type: TabStopType.RIGHT,  position: CONTENT_W },
    ],
    children: [
      new TextRun({ text: '\t', font: TIMES, size: SIZE_BODY }),
      mathNode,
      new TextRun({ text: '\t' + label, font: TIMES, size: SIZE_BODY }),
    ],
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
// Inline-формула (без centered display-режима) — для встраивания внутрь
// абзаца с TextRun'ами.
function inlineMath(...nodes) {
  return new MathBlock({ children: nodes.map(_m) });
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
// Bookmark ID для главы, чтобы в оглавлении делать на неё кликабельную ссылку.
// num: '1'..'11', 'А', 'Б' (для приложений), '' (для безномерных вроде «Введение»).
function bookmarkIdFor(num, fallback) {
  if (!num) return 'sec_' + (fallback ? String(fallback).toLowerCase() : 'unnumbered');
  // Приложения по-русски: А, Б -> a, b
  const map = { 'А': 'a', 'Б': 'b' };
  return 'sec_' + (map[num] || num);
}
// Заголовок раздела (1, 2, ...) — прописными, по центру, жирный, по ГОСТ 7.32
function h1(num, text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    alignment: AlignmentType.CENTER,
    spacing: { before: 480, after: 240, line: 360, lineRule: 'auto' },
    pageBreakBefore: true,
    children: [
      new Bookmark({
        id: bookmarkIdFor(num),
        children: [
          new TextRun({ text: `${num} ${text}`.toUpperCase(),
                        font: TIMES, size: SIZE_BODY, bold: true }),
        ],
      }),
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
// Заголовок без номера (Содержание, Введение, Приложение А) — ПРОПИСНЫМИ, по центру.
// opts.id — bookmark ID, на который будет ссылаться оглавление.
function h1Unnumbered(text, opts = {}) {
  const id = opts.id ? bookmarkIdFor('', opts.id) : null;
  const titleRun = new TextRun({ text: text.toUpperCase(), font: TIMES, size: SIZE_BODY, bold: true });
  const titleNode = id
    ? new Bookmark({ id, children: [titleRun] })
    : titleRun;
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    alignment: AlignmentType.CENTER,
    spacing: { before: 480, after: 360, line: 360, lineRule: 'auto' },
    pageBreakBefore: opts.pageBreakBefore !== false,
    children: [ titleNode ],
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
// Маркированный пункт с математическим выражением в начале:
//   mdash(varExpr, prose)
//   mdash([sub('γ', 'n')], 'Лоренц-фактор на шаге n.')
// Первый аргумент — строка (станет MathRun), один Math-узел или массив узлов.
// Второй — обычный поясняющий текст, идущий после « — ».
function mdash(varExpr, prose) {
  let nodes;
  if (typeof varExpr === 'string')      nodes = [_m(varExpr)];
  else if (Array.isArray(varExpr))      nodes = varExpr.map(_m);
  else                                  nodes = [varExpr];
  return new Paragraph({
    spacing: { line: 360, lineRule: 'auto', after: 60 },
    indent: { left: mm(7.5), hanging: mm(5) },
    alignment: AlignmentType.JUSTIFIED,
    children: [
      new TextRun({ text: '— ', font: TIMES, size: SIZE_BODY }),
      new MathBlock({ children: nodes }),
      new TextRun({ text: ' — ' + prose, font: TIMES, size: SIZE_BODY }),
    ],
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
// Кликабельный пункт оглавления. anchorId — bookmark в теле документа,
// на который ссылка должна указывать. Номер страницы оставлен пустым —
// проставляется вручную в Word после генерации.
function tocItem(anchorId, num, text) {
  const tabs = [{ type: TabStopType.RIGHT, position: CONTENT_W - mm(2),
                  leader: 'dot' }];
  return new Paragraph({
    spacing: { line: 360, lineRule: 'auto', after: 60 },
    indent: { firstLine: 0 },
    tabStops: tabs,
    children: [
      new InternalHyperlink({
        anchor: anchorId,
        children: [
          new TextRun({ text: `${num ? num + ' ' : ''}${text}`,
                        font: TIMES, size: SIZE_BODY }),
          new TextRun({ text: '\t', font: TIMES, size: SIZE_BODY }),
        ],
      }),
    ],
  });
}

// Состав документа (для содержания): [anchorId, num, text].
// anchorId должен совпадать с bookmark'ом, который ставит h1 / h1Unnumbered.
const TOC_ENTRIES = [
  ['sec_intro', '',   'Введение'],
  ['sec_1',  '1',  'Интегратор Бориса — релятивистский симплектический интегратор'],
  ['sec_2',  '2',  'Трилинейная интерполяция поля с учётом YZ-симметрии'],
  ['sec_3',  '3',  'Двумерный bucket-индекс магнитов'],
  ['sec_4',  '4',  'Вклад магнита с прескринингом по охватывающей сфере'],
  ['sec_5',  '5',  'Детектор пересечения финишной плоскости'],
  ['sec_6',  '6',  'Расстояние от точки до полилинии проектной орбиты'],
  ['sec_7',  '7',  'Составная функция потерь оптимизатора'],
  ['sec_8',  '8',  'Стохастический оптимизатор (hill-climber c киками)'],
  ['sec_9',  '9',  'Генерация решётки из LTE-описания'],
  ['sec_10', '10', 'Сборка таблицы конфигурации и суб-решёток'],
  ['sec_11', '11', 'Анализ огибающей и бетатронных частот по FFT'],
  ['sec_a',  'А',  'Приложение А. Файлы конфигурации'],
  ['sec_b',  'Б',  'Приложение Б. Файлы данных'],
];

// ---------------------------------------------------------------------------
// Введение
// ---------------------------------------------------------------------------
startSection(null);
const intro = [
  h1Unnumbered('Введение', { id: 'intro', pageBreakBefore: true }),
  body('Настоящий документ содержит описание вычислительных моделей, реализованных в составе программ SynchrotronTracker.exe и SynchrotronOptimizer.exe. Документ оформлен в соответствии с методическими указаниями РД 50-34.698-90, разделом «Описание алгоритма», где для каждого алгоритма приводятся пять подразделов: назначение и характеристика, используемая информация, результаты решения, математическое описание, алгоритм решения.'),
  body('Источниками исходных данных служат файлы решётки в формате Elegant LTE (workspace/lte/), бинарные карты поля Ansys (workspace/fields_bin/), JSON-конфигурации (workspace/configs/) и CSV-таблицы базы данных (workspace/tables/). Результатами работы программ являются бинарные файлы траекторий (workspace/results/trajectories_exp{id}.bin) и журналы оптимизации в формате JSON.'),
  body('Все размерные величины приведены в системе СИ. Ссылки на исходный код имеют вид «файл:строка», например, src/Integrator.h:53. Документ описывает 11 алгоритмов и содержит два приложения со справочной информацией о конфигурационных файлах и форматах данных.'),
];

// ---------------------------------------------------------------------------
// 1. Boris pusher
// ---------------------------------------------------------------------------
startSection('1');
const alg1 = [
  h1('1', 'Интегратор Бориса — релятивистский симплектический интегратор'),
  body('Исходный код: src/Integrator.h:53, метод Integrator::borisStep.'),

  h2('1.1', 'Назначение и характеристика'),
  body('Интегратор Бориса представляет собой одношаговый релятивистский метод продвижения заряженной частицы в произвольном электромагнитном поле. В данной работе применяется его реализация по схеме расщепления drift–kick–drift. Метод обладает следующими свойствами:'),
  dash('симплектичен — сохраняет фазовый объём;'),
  dash('точно сохраняет модуль импульса |p|² для чисто магнитного поля, поскольку шаг поворота является ортогональным преобразованием по построению;'),
  dash('симметричен по времени, формальный порядок аппроксимации составляет O(dt²);'),
  dash('требует одного вычисления поля за шаг — против четырёх в методе Рунге–Кутты 4-го порядка, что даёт ускорение приблизительно в 4 раза при идентичном бюджете шагов.'),
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
  body('Движение релятивистской заряженной частицы в электромагнитном поле описывается классическими уравнениями Лоренца:'),
  mathPara(frac('dx', 'dt'), ' = v,    ', frac('dp', 'dt'), ' = q (E + v × B)'),
  body('Здесь приняты следующие обозначения:'),
  mdash('x(t)', 'координаты частицы в глобальной системе координат, м.'),
  mdash('v(t)', 'скорость частицы в глобальной системе координат, м/с.'),
  mdash('p(t)', 'релятивистский импульс, кг·м/с — величина, непосредственно интегрируемая в коде; совпадает с классическим mv только при малых скоростях.'),
  mdash('q', 'заряд частицы, Кл; для электронов принимается q = −1,602·10⁻¹⁹ Кл.'),
  mdash('m', 'масса покоя, кг; для электронов принимается m = 9,109·10⁻³¹ кг.'),
  mdash(['E,  B'], 'векторы электрического и магнитного полей в точке частицы, В/м и Тл соответственно.'),
  mdash('t', 'физическое (лабораторное) время, с.'),
  body('Стоит отметить, что векторное произведение v × B соответствует магнитной составляющей силы Лоренца, отнесённой к единице заряда. Магнитное поле искривляет траекторию частицы, однако работы не совершает, и кинетическая энергия сохраняется точно.'),
  body('Связь импульса и скорости задаётся через релятивистский Лоренц-фактор:'),
  mathPara('p = γ m v'),
  mathPara(
    'γ = ',
    frac('1', rad(['1 − ', frac(sup('|v|', '2'), sup('c', '2'))])),
    ' = ',
    rad(['1 + ', frac(sup('|p|', '2'), [sup('m', '2'), ' ', sup('c', '2')])]),
  ),
  body('где:'),
  mdash('γ ≥ 1', 'Лоренц-фактор, безразмерная величина, характеризующая степень релятивистичности движения (γ = 1 в покое, γ → ∞ при |v| → c).'),
  mdash('c = 299 792 458 м/с', 'скорость света в вакууме (точное значение по определению метра).'),
  mdash([sup('|v|', '2'), ' = ', sup(sub('v', 'x'), '2'), ' + ', sup(sub('v', 'y'), '2'), ' + ', sup(sub('v', 'z'), '2')], 'квадрат модуля скорости.'),
  mdash([sup('|p|', '2'), ' = ', sup(sub('p', 'x'), '2'), ' + ', sup(sub('p', 'y'), '2'), ' + ', sup(sub('p', 'z'), '2')], 'квадрат модуля импульса.'),
  body('Стоит отметить, что в коде применяется именно вторая форма — через |p|: значение γ восстанавливается непосредственно по импульсу, и хранить скорость отдельно не требуется. В рассматриваемой решётке все устройства являются магнитными, поэтому E ≡ 0, и шаги E-толчка, рассмотренные ниже, не оказывают воздействия на импульс. Тем не менее общая форма уравнений сохранена — для возможного включения электростатических элементов в будущем.'),
  body('Алгоритм Бориса представляет собой расщепление одного шага Δt на семь последовательных подопераций над состоянием частицы (xₙ, pₙ) → (xₙ₊₁, pₙ₊₁), которые описаны далее.'),
  body('На первом этапе выполняется полудрейф — координаты частицы продвигаются на половину шага по текущему импульсу:'),
  mathPara(sub('x', 'n+½'), ' = ', sub('x', 'n'), ' + ',
           frac(sub('p', 'n'), [sub('γ', 'n'), ' m']),
           ' · ', frac('Δt', '2')),
  body('где:'),
  mdash('Δt', 'шаг интегрирования, с; задаётся в столбце dt таблицы experiments.csv.'),
  mdash([sub('γ', 'n'), ' = ', rad(['1 + ', frac(sup(sub('|p|', 'n'), '2'), [sup('m', '2'), ' ', sup('c', '2')])])], 'Лоренц-фактор на шаге n.'),
  body('На втором этапе производится опрос суммарного электромагнитного поля в полученной середине шага:'),
  mathPara('(B, E) ← lattice.getTotalFields(', sub('x', 'n+½'), ')'),
  body('Следует подчеркнуть, что данный вызов является единственным обращением к карте поля за один шаг интегрирования; метод Рунге–Кутты 4-го порядка потребовал бы четырёх таких вычислений.'),
  body('На третьем этапе выполняется первый полу-E-толчок:'),
  mathPara(sup('p', '−'), ' = ', sub('p', 'n'), ' + ', frac('q E Δt', '2')),
  body('где:'),
  mdash(sup('p', '−'), 'промежуточный импульс после первого полу-E-толчка. При E = 0 справедливо p⁻ = pₙ.'),
  body('На четвёртом этапе пересчитывается Лоренц-фактор по полученному импульсу:'),
  mathPara(sup('γ', '−'), ' = ', rad(['1 + ', frac(sup(['|', sup('p', '−'), '|'], '2'), [sup('m', '2'), ' ', sup('c', '2')])])),
  body('Полученное γ⁻ представляет собой мгновенный релятивистский фактор, используемый на последующем магнитном повороте.'),
  body('На пятом этапе выполняется собственно магнитный поворот импульса. Для этого предварительно формируются два вспомогательных вектора:'),
  mathPara('t = ', frac('q B Δt', ['2 ', sup('γ', '−'), ' m'])),
  mathPara('s = ', frac('2 t', ['1 + ', sup('|t|', '2')])),
  body('где:'),
  mdash('t', 'половинный вектор поворота; его модуль равен tan(θ/2), где θ = qB·Δt/(γ⁻·m) — угол ларморовского поворота за один шаг.'),
  mdash('s', 'перенормированный вектор поворота, подобранный таким образом, чтобы двухшаговое обновление импульса, рассматриваемое далее, точно сохраняло модуль |p|.'),
  body('Сам поворот импульса производится в два этапа:'),
  mathPara('p′ = ', sup('p', '−'), ' + ', sup('p', '−'), ' × t'),
  mathPara(sup('p', '+'), ' = ', sup('p', '−'), ' + p′ × s,    |', sup('p', '+'), '| = |', sup('p', '−'), '| точно'),
  body('где:'),
  mdash('p′', 'промежуточный импульс, полученный после первого полуповорота.'),
  mdash(sup('p', '+'), 'импульс после полного магнитного поворота. Равенство |p⁺| = |p⁻| справедливо как алгебраическое тождество — именно благодаря этому свойству интегратор Бориса точно сохраняет энергию частицы в чисто магнитном поле.'),
  body('На шестом этапе выполняется второй полу-E-толчок:'),
  mathPara(sub('p', 'n+1'), ' = ', sup('p', '+'), ' + ', frac('q E Δt', '2')),
  body('Совместно с третьим этапом этот толчок реализует центрированное E-ускорение: каждый из полутолчков вносит в импульс величину q·E·Δt/2.'),
  body('Наконец, на седьмом этапе производится второй полудрейф, завершающий шаг интегрирования:'),
  mathPara(sub('x', 'n+1'), ' = ', sub('x', 'n+½'), ' + ',
           frac(sub('p', 'n+1'), [sub('γ', 'n+1'), ' m']),
           ' · ', frac('Δt', '2')),
  body('где:'),
  mdash([sub('γ', 'n+1'), ' = ', rad(['1 + ', frac(sup(sub('|p|', 'n+1'), '2'), [sup('m', '2'), ' ', sup('c', '2')])])], 'Лоренц-фактор на новом шаге.'),
  body('Таким образом, по итогам выполнения семи подопераций состояние частицы продвинуто ровно на один шаг Δt с формальным порядком аппроксимации O(Δt²).'),

  h2('1.5', 'Алгоритм решения'),
  code(
    'Вход:  s = (pos, mom), m0, q, lattice, dt\n' +
    'Выход: s обновлено до шага n+1\n\n' +
    ' 1.  γ ← sqrt(1 + |mom|² / (m0²·c²))            ; Лоренц-фактор по текущему импульсу\n' +
    ' 2.  pos ← pos + mom · (dt/2) / (m0·γ)          ; ПЕРВЫЙ ПОЛУДРЕЙФ: x_n → x_{n+½}\n' +
    ' 3.  (B, E) ← lattice.getTotalFields(pos)       ; ОПРОС ПОЛЯ в середине шага\n' +
    '                                                ;   (единственное обращение к карте за шаг)\n' +
    ' 4.  p_minus ← mom + E · q · dt/2               ; ПЕРВЫЙ ПОЛУ-E-ТОЛЧОК (при E=0 — no-op)\n' +
    ' 5.  γ_minus ← sqrt(1 + |p_minus|² / (m0²·c²))  ; пересчёт γ по промежуточному импульсу\n' +
    ' 6.  t ← B · q · dt / (2·m0·γ_minus)            ; полуугольный вектор поворота\n' +
    ' 7.  s ← 2·t / (1 + |t|²)                       ; перенормированный вектор поворота\n' +
    ' 8.  p_prime ← p_minus + cross(p_minus, t)      ; первый полуповорот\n' +
    ' 9.  p_plus  ← p_minus + cross(p_prime, s)      ; завершение поворота; |p_plus| = |p_minus|\n' +
    '10.  mom ← p_plus + E · q · dt/2                ; ВТОРОЙ ПОЛУ-E-ТОЛЧОК\n' +
    '11.  γ ← sqrt(1 + |mom|² / (m0²·c²))            ; γ на новом шаге\n' +
    '12.  pos ← pos + mom · (dt/2) / (m0·γ)          ; ВТОРОЙ ПОЛУДРЕЙФ: x_{n+½} → x_{n+1}'),
  body('В ранних версиях существовали вспомогательные методы Integrator::driftStep и AcceleratorConfig::isDriftRegion для замены полного шага Бориса аналитическим дрейфом x += v·dt на участках без магнитов. Оба удалены в рефакторинге 2026-06: прескрининг по bucket-индексу (см. раздел 3) уже сводит шаг в пустой области к одному обращению к бакету и пустому циклу по кандидатам, и отдельная ветка для дрейфа практической пользы не давала.'),
  body('В текущей конфигурации все устройства магнитные, поэтому E ≡ 0 и шаги 3, 6 (полу-E-толчки) не изменяют импульс. Алгоритм оставлен в общем виде для совместимости при возможном появлении электрических элементов.'),
];

// ---------------------------------------------------------------------------
// 2. Trilinear field interpolation
// ---------------------------------------------------------------------------
startSection('2');
const alg2 = [
  h1('2', 'Трилинейная интерполяция поля с учётом YZ-симметрии'),
  body('Исходный код: src/FieldMap.h:41, метод FieldMap::getInterpolatedField.'),

  h2('2.1', 'Назначение и характеристика'),
  body('Метод трилинейной интерполяции применяется для вычисления вектора магнитного (или электрического) поля в произвольной точке локальной системы координат магнита по предварительно рассчитанной карте поля. Карта представляет собой регулярную декартову сетку, полученную как результат моделирования в Ansys. Трилинейная интерполяция обеспечивает непрерывность C⁰ и совпадение со значением в узле сетки. Стоит отметить, что опциональное зеркалирование по плоскости YZ (флаг useSymmetry) позволяет хранить только половину объёма для полей, антисимметричных по координате x, что вдвое снижает требования к памяти.'),

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
  body('Возвращаемое значение представляет собой интерполированный вектор поля FieldVector в локальных координатах. Для точек, лежащих за пределами ограничивающего параллелепипеда сетки, возвращается нулевой вектор (0, 0, 0).'),

  h2('2.4', 'Математическое описание'),
  body('Карта поля задана в узлах регулярной декартовой сетки. Точка опроса p = (px, py, pz), заданная в локальной системе координат, попадает в одну из единичных ячеек сетки; вектор поля в этой точке определяется как взвешенное среднее значений в восьми углах ячейки.'),
  body('Дробные координаты ячейки рассчитываются по формулам:'),
  mathPara('i = ⌊', frac([sub('p', 'x'), ' − ', sub('x', 'min')], 'Δx'), '⌋,    ',
           sub('t', 'x'), ' = ',
           frac([sub('p', 'x'), ' − ', sub('x', 'min')], 'Δx'), ' − i'),
  mathPara('j = ⌊', frac([sub('p', 'y'), ' − ', sub('y', 'min')], 'Δy'), '⌋,    ',
           sub('t', 'y'), ' = ',
           frac([sub('p', 'y'), ' − ', sub('y', 'min')], 'Δy'), ' − j'),
  mathPara('k = ⌊', frac([sub('p', 'z'), ' − ', sub('z', 'min')], 'Δz'), '⌋,    ',
           sub('t', 'z'), ' = ',
           frac([sub('p', 'z'), ' − ', sub('z', 'min')], 'Δz'), ' − k'),
  body('где:'),
  mdash(['i,  j,  k'], 'целые индексы нижнего угла ячейки (нумерация с нуля; верхний угол ячейки имеет индексы (i+1, j+1, k+1)).'),
  mdash([sub('t', 'x'), ',  ', sub('t', 'y'), ',  ', sub('t', 'z'), ' ∈ [0, 1]'], 'дробные смещения внутри ячейки; значение t_α = 0 соответствует нижней грани, t_α = 1 — верхней.'),
  mdash([sub('x', 'min'), ',  ', sub('y', 'min'), ',  ', sub('z', 'min')], 'координаты угла сетки с наименьшими индексами, м.'),
  mdash(['Δx,  Δy,  Δz'], 'шаг сетки по соответствующим осям, м (равномерный в пределах одной оси, но независимый между осями).'),
  body('Введём обозначение v_{abc} для вектора поля в угле с индексами (i+a, j+b, k+c), где a, b, c ∈ {0, 1}. Тогда интерполированное значение задаётся выражением:'),
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
  body('где:'),
  mdash(sub('v', 'abc'), 'значение поля в одном из восьми углов ячейки (Тл — для магнитного поля B, В/м — для электрического E).'),
  mdash(sub('w', 'abc'), 'неотрицательный весовой коэффициент, линейно убывающий с расстоянием от точки опроса до соответствующего угла. Сумма всех восьми весов равна единице, и результат, таким образом, представляет собой истинное взвешенное среднее.'),
  body('Следует отметить, что в угловых точках (t_α ∈ {0, 1} по всем осям) ровно один весовой коэффициент равен единице, а остальные обращаются в ноль. Следовательно, F принимает в узлах сетки точные табличные значения — интерполянт согласован в узлах.'),
  body('Правило YZ-симметрии применяется в том случае, когда карта поля антисимметрична по координате x и хранится только половина объёма (x ∈ [0, x_max]). Запрос с отрицательным x сворачивается в сохранённую половину по правилу:'),
  mathPara('F(', sub('p','x'),', ', sub('p','y'),', ', sub('p','z'), ') = ',
           sub('F', 'raw'), '(', sub('p','x'),', ', sub('p','y'),', ', sub('p','z'), '),    если ', sub('p','x'), ' ≥ 0'),
  mathPara('F(', sub('p','x'),', ', sub('p','y'),', ', sub('p','z'), ') = ',
           'diag(−1, +1, +1) · ', sub('F', 'raw'), '(−', sub('p','x'),', ', sub('p','y'),', ', sub('p','z'), '),    если ', sub('p','x'), ' < 0'),
  body('где:'),
  mdash(sub('F', 'raw'), 'результат интерполяции, выполненной на сохранённой (px ≥ 0) половине сетки.'),
  mdash('diag(−1, +1, +1)', 'диагональная матрица, инвертирующая x-компоненту и оставляющая остальные без изменений. Применение данной матрицы отражает физическую антисимметрию компоненты B_x по плоскости YZ — тогда как компоненты B_y и B_z остаются симметричными.'),

  h2('2.5', 'Алгоритм решения'),
  code(
    'Вход:  p, useSymmetry, grid (header + data)\n' +
    'Выход: вектор поля в точке p\n\n' +
    '1. mirrored ← false                            ; флаг применения YZ-зеркалирования\n' +
    '2. если useSymmetry и p.x < 0:                 ; сворачиваем запрос в сохранённую половину сетки\n' +
    '       p.x ← -p.x; mirrored ← true\n' +
    '3. если p вне [minX..maxZ]: вернуть {0,0,0}    ; точка за пределами bbox → поле = 0\n' +
    '4. fx ← (p.x - minX)/stepX;  i ← floor(fx);  tx ← fx - i      ; дробные координаты ячейки\n' +
    '   fy ← (p.y - minY)/stepY;  j ← floor(fy);  ty ← fy - j\n' +
    '   fz ← (p.z - minZ)/stepZ;  k ← floor(fz);  tz ← fz - k\n' +
    '5. ограничить i, j, k интервалом [0, n-2]      ; защита для граничного случая p == max\n' +
    '6. извлечь 8 значений v000..v111 из data[]     ; угловые значения текущей ячейки\n' +
    '7. покомпонентно вычислить F по формуле трилинейной интерполяции (см. § 2.4)\n' +
    '8. если mirrored: F.x ← -F.x                   ; инверсия x-компоненты для YZ-симметрии\n' +
    '9. вернуть F'),
];

// ---------------------------------------------------------------------------
// 3. Spatial bucket grid
// ---------------------------------------------------------------------------
startSection('3');
const alg3 = [
  h1('3', 'Двумерный bucket-индекс магнитов'),
  body('Исходный код: src/LatticeManager.h:179 (buildSpatialIndex), src/LatticeManager.h:235 (getTotalFields).'),

  h2('3.1', 'Назначение и характеристика'),
  body('Двумерный bucket-индекс применяется для ускорения поиска магнитов, потенциально вносящих вклад в поле в точке текущей позиции частицы. Без подобного индекса симулятор был бы вынужден перебирать все ~2762 элемента кольца на каждом шаге интегрирования. Применение индекса позволяет сократить число кандидатов до ~5 на шаг; на участках чистого дрейфа соответствующий бакет оказывается пустым, и стоимость шага сводится к единственному обращению к индексу без каких-либо вычислений поля.'),
  body('Индекс строится однократно — после завершения добавления элементов через метод addElement(...) — и далее используется в режиме «только чтение».'),

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
  body('Каждому магниту сопоставлена заранее вычисленная ограничивающая сфера с центром c_i (в глобальных координатах) и радиусом r_i, полностью покрывающая его карту поля. Индекс строится в плоскости XY, поскольку кольцо синхротрона лежит в плоскости z ≈ 0.'),
  body('Размер бакета определяется по максимальному радиусу ограничивающей сферы среди всех элементов. Шаг сетки задаётся выражением:'),
  mathPara('h = max(1 м,  2 ', sub('r', 'max'), ')'),
  body('где:'),
  mdash('h', 'длина ребра одного квадратного бакета, м.'),
  mdash([sub('r', 'max'), ' = ', sub('max', 'i'), ' ', sub('r', 'i')], 'максимальный радиус ограничивающей сферы по всем магнитам, м. Множитель 2 гарантирует, что XY-проекция любой ограничивающей окружности укладывается в группу 2×2 бакетов, и магнит регистрируется не более чем в четырёх бакетах.'),
  dash('Нижняя граница в 1 метр предотвращает образование патологически мелких бакетов при наборе мелких магнитов — в противном случае массив бакетов разрастался бы непропорционально.'),
  body('Каждой точке (x, y, z) в глобальной системе координат сопоставляются индексы бакета (b_x, b_y) по формулам:'),
  mathPara(sub('b', 'x'), ' = ⌊', frac('x', 'h'), '⌋ − ', subSup('b', 'x', 'min')),
  mathPara(sub('b', 'y'), ' = ⌊', frac('y', 'h'), '⌋ − ', subSup('b', 'y', 'min')),
  body('где:'),
  mdash('⌊·⌋', 'операция округления вниз до ближайшего целого (floor).'),
  mdash([subSup('b', 'x', 'min'), ',  ', subSup('b', 'y', 'min')], 'глобальные целочисленные смещения, подобранные на этапе построения индекса таким образом, чтобы все занятые бакеты имели неотрицательные индексы. Значения вычисляются по наименьшим x и y, достигаемым ограничивающими сферами.'),
  body('Таким образом, магнит i регистрируется в каждом бакете (b_x, b_y), пересекающем диск { (x, y) : (x − c_{i,x})² + (y − c_{i,y})² ≤ r_i² }. По построению число таких бакетов не превышает 2×2 = 4.'),

  h2('3.5', 'Алгоритм решения'),
  body('Построение индекса (выполняется однократно):'),
  code(
    '1. maxR ← max_i elements[i].worldRadius        ; наибольший радиус охватывающей сферы\n' +
    '2. bucketSize ← max(1.0, 2 · maxR)             ; шаг сетки = 2·r_max, но не менее 1 м\n' +
    '                                                ; нижняя граница защищает от мелких бакетов\n' +
    '3. xLo, xHi ← мин/макс по элементам от (worldCenter.x ± worldRadius)\n' +
    '   yLo, yHi ← аналогично по y                   ; bbox занимаемой кольцом области в XY\n' +
    '4. bxMin ← floor(xLo/bucketSize);  byMin ← floor(yLo/bucketSize)   ; смещения индексов\n' +
    '   bxN ← floor(xHi/bucketSize) - bxMin + 1                          ; размеры сетки\n' +
    '   byN ← floor(yHi/bucketSize) - byMin + 1\n' +
    '5. bucketGrid ← массив пустых списков размера bxN·byN              ; плоский 2D-массив\n' +
    '6. для каждого элемента i:                     ; регистрируем магниты в бакетах\n' +
    '       bxLo..bxHi ← диапазон бакетов по x      ; XY-проекция охватывающей сферы\n' +
    '       byLo..byHi ← диапазон бакетов по y      ; (по построению не более 2×2 бакетов)\n' +
    '       для bx ∈ [bxLo..bxHi], by ∈ [byLo..byHi]:\n' +
    '           bucketGrid[bx·byN + by].push_back(i)'),
  body('Запрос (на каждом шаге частицы):'),
  code(
    '1. bx ← floor(pos.x/bucketSize) - bxMin        ; индексы бакета, в который попала частица\n' +
    '   by ← floor(pos.y/bucketSize) - byMin\n' +
    '2. если (bx, by) вне сетки: B, E ← 0; вернуть  ; вышли за пределы кольца → поля нет\n' +
    '3. для каждого индекса i из bucketGrid[bx·byN + by]:   ; перебор только локальных кандидатов\n' +
    '       elements[i].addContribution(pos, B, E)          ; вклад магнита (с прескринингом)'),
  body('На участках дрейфа бакет пуст, поэтому внутренний цикл по кандидатам не выполняется и вычислений поля не происходит.'),
];

// ---------------------------------------------------------------------------
// 4. Magnet contribution
// ---------------------------------------------------------------------------
startSection('4');
const alg4 = [
  h1('4', 'Вклад магнита с прескринингом по охватывающей сфере'),
  body('Исходный код: src/LatticeManager.h:114, метод MagneticElement::addContribution.'),

  h2('4.1', 'Назначение и характеристика'),
  body('Метод addContribution используется для добавления вклада одного магнита в накапливающиеся суммы полей (B, E) в точке опроса, заданной в глобальной системе координат. Перед выполнением полного преобразования координат и интерполяции производится быстрая проверка по охватывающей сфере (7 арифметических операций), позволяющая отсечь заведомо далёкие элементы и тем самым существенно снизить вычислительные затраты.'),

  h2('4.2', 'Используемая информация'),
  tableCaption('4.1', 'Входные данные алгоритма'),
  table3col(
    ['Обозначение', 'Смысл', 'Источник'],
    [
      ['globalPos', 'точка опроса в глобальных координатах', 'цикл трекинга'],
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
  body('Метод addContribution маршрутизирует вклад в зависимости от типа поля: значения с типом MAGNETIC накапливаются в накопителе totalB, ELECTRIC — в totalE. Стоит отметить, что все текущие устройства являются магнитными, поэтому totalE на практике остаётся равным нулю, а соответствующие E-толчки в интеграторе Бориса не оказывают воздействия. Тем не менее логика поддерживает оба типа поля, и добавление электростатического элемента в будущем не потребует изменения кода.'),

  h2('4.3', 'Результаты решения'),
  body('Результатом работы метода является обновление по месту накопителя totalB (если type == MAGNETIC) или totalE (если type == ELECTRIC). Возвращаемое значение не предусмотрено.'),

  h2('4.4', 'Математическое описание'),
  body('Каждому магниту сопоставлена собственная локальная система координат, в которой и вычислена его карта поля. Вычисление вклада магнита i в произвольной точке x глобальной системы координат складывается из четырёх геометрических шагов и заключительного масштабирования.'),
  body('На первом этапе выполняется прескрининг по охватывающей сфере. Магнит сразу отбрасывается, если точка опроса находится за пределами его сферы влияния, что определяется условием:'),
  mathPara(sup(['|x − ', sub('c', 'i'), '|'], '2'), ' > ', sup(sub('r', 'i'), '2'),
           '    ⇒    вклад = 0'),
  body('где:'),
  mdash(sub('c', 'i'), 'центр охватывающей сферы магнита i в глобальных координатах, м.'),
  mdash(sub('r', 'i'), 'радиус охватывающей сферы, м.'),
  body('Следует подчеркнуть, что возведение в квадрат с обеих сторон неравенства позволяет обойтись без вычисления корня в горячем цикле трекинга.'),
  body('На втором этапе выполняется преобразование точки опроса из глобальной системы координат в локальную систему магнита:'),
  mathPara(sub('x', 'loc'), ' = ', sup(sub('R', 'i'), '⊤'), ' (x − ', sub('T', 'i'), ')'),
  body('где:'),
  mdash(sub('T', 'i'), 'трансляция, задающая положение локального начала координат магнита в глобальной системе, м.'),
  mdash(sub('R', 'i'), 'матрица поворота 3×3 (ортогональная, удовлетворяющая RᵀR = I), переводящая локальный вектор в глобальную систему. Матрица строится по углам yaw, pitch, roll из таблицы lattice_configs.csv.'),
  mdash('Rᵀ', 'обратное преобразование. Поскольку для ортогональной матрицы транспонирование совпадает с обращением, операция возвращает точку из глобальной системы в локальную.'),
  body('На третьем этапе по полученной локальной координате производится интерполяция поля по карте, выполненная по методу § 2.4:'),
  mathPara(sub('F', 'loc'), ' = trilinear(', sub('map', 'i'), ', ', sub('x', 'loc'), ', ', sub('useSymmetry', 'i'), ')'),
  body('где:'),
  mdash(sub('map', 'i'), 'карта поля магнита i, представленная в виде двоичной сетки векторов.'),
  mdash([sub('useSymmetry', 'i'), ' ∈ {true, false}'], 'признак применения YZ-зеркалирования (см. § 2.4).'),
  body('В случае, если F_loc = 0 (точка опроса попадает в нулевую зону сетки), дальнейшие вычисления прекращаются и метод сразу завершает работу.'),
  body('На четвёртом этапе выполняется поворот вектора поля обратно в глобальную систему координат:'),
  mathPara(sub('F', 'world'), ' = ', sub('R', 'i'), ' · ', sub('F', 'loc')),
  body('Трансляция на данном этапе не применяется. Магнитное поле представляет собой вектор, а не точку пространства, поэтому единственной действующей частью преобразования является поворот.'),
  body('Наконец, на пятом этапе производится масштабирование. Все индивидуальные калибровочные множители магнита перемножаются по формуле:'),
  mathPara(sub('scale', 'i'), ' = ',
           subSup('c', 'i', 'hw'), ' · ', subSup('c', 'i', 'sp'), ' · ', sub('d', 'i'),
           ' · '),
  mathPara('     · ',
           '{ ', frac(sub('I', 'i'), subSup('I', 'i', 'nom')), '   если магнит i активен;    1   иначе }'),
  body('где:'),
  mdash(subSup('c', 'i', 'hw'), 'коэффициент калибровки на уровне устройства (hardware_library.correct_coeff), безразмерный; задаётся по одному значению на тип устройства.'),
  mdash(subSup('c', 'i', 'sp'), 'индивидуальный коэффициент калибровки (sp_elements.correct_coeff), безразмерный; задаётся по одному значению на пару (sp_type, sp_element_id), по умолчанию принимается равным 1.0.'),
  mdash([sub('d', 'i'), ' ∈ {+1, −1}'], 'признак инверсии направления (hardware_library.field_direction). Применяется в тех случаях, когда физическая ориентация магнита отличается от заложенной в карту поля.'),
  mdash(sub('I', 'i'), 'фактический ток через магнит, А (sp_elements.arg_val). Используется только для активных электромагнитов; для постоянных магнитов тока нет и множитель принимается равным 1.'),
  mdash(subSup('I', 'i', 'nom'), 'номинальный ток, при котором рассчитывалась карта поля, А (hardware_library.nominal_arg). Отношение I_i / I_i^nom линейно масштабирует поле по току — в предположении ненасыщенного режима.'),
  body('На шестом, заключительном этапе накопленное значение добавляется к соответствующей сумме поля:'),
  mathPara(sub('B', 'total'), ' += ', sub('scale', 'i'), ' · ', sub('F', 'world'), ',    если ', sub('type', 'i'), ' = MAGNETIC'),
  mathPara(sub('E', 'total'), ' += ', sub('scale', 'i'), ' · ', sub('F', 'world'), ',    если ', sub('type', 'i'), ' = ELECTRIC'),
  body('где:'),
  mdash(['totalB,  totalE'], 'текущие суммы магнитного и электрического полей в точке опроса, Тл и В/м соответственно.'),
  mdash(sub('type', 'i'), 'тип поля: MAGNETIC (по умолчанию для всех текущих устройств) или ELECTRIC (заготовка для будущих электростатических элементов). Значение определяется по колонке field_type таблицы hardware_library.csv.'),

  h2('4.5', 'Алгоритм решения'),
  code(
    '1. dx,dy,dz ← globalPos − worldCenter           ; ПРЕСКРИНИНГ по охватывающей сфере:\n' +
    '   если dx²+dy²+dz² > worldRadius²:  return     ; точка вне сферы влияния → вклад нулевой\n' +
    '                                                ; (сравнение квадратов, без sqrt)\n' +
    '2. localPos ← transform.transformPointInv(globalPos)   ; ГЛОБАЛЬНАЯ → ЛОКАЛЬНАЯ система\n' +
    '3. localF   ← map.getInterpolatedField(localPos, useSymmetry)\n' +
    '                                                ; интерполяция по карте поля (см. § 2)\n' +
    '4. если localF == 0:  return                    ; нулевая зона карты — нет смысла продолжать\n' +
    '5. globalF ← transform.rotateVector(localF)     ; ЛОКАЛЬНАЯ → ГЛОБАЛЬНАЯ (только поворот)\n' +
    '6. scale ← correctCoeff · fieldDirection        ; калибровка по типу + знак направления\n' +
    '   если isActive: scale *= currentScale         ; масштаб по току для активных магнитов\n' +
    '7. если type == MAGNETIC: totalB += scale · globalF    ; маршрутизация по типу поля\n' +
    '   иначе                  totalE += scale · globalF'),
];

// ---------------------------------------------------------------------------
// 5. Finish plane
// ---------------------------------------------------------------------------
startSection('5');
const alg5 = [
  h1('5', 'Детектор пересечения финишной плоскости'),
  body('Исходный код: src/FinishPlane.h:51 (FinishPlane::hasCrossed), setup на строке 30.'),

  h2('5.1', 'Назначение и характеристика'),
  body('Детектор используется для определения момента, когда частица пересекает ориентированную плоскость, заданную тремя точками. Каждая частица обладает собственным состоянием детектора — счётчиком вызовов и последним знаковым расстоянием, — что обеспечивает корректную работу в многопоточном режиме OpenMP. Стоит отметить, что параметр skipSteps позволяет пропускать заданное число первых вызовов. Данный механизм необходим в тех случаях, когда плоскость пересекает орбиту несколько раз: например, плоскость, поставленная на половине кольца, проходит и вблизи стартовой точки, и около цели.'),

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
  body('Метод hasCrossed(pos) возвращает значение true при первом вызове, в котором знаковое расстояние от точки pos до плоскости сменило знак относительно сохранённого значения lastD. Полученный флаг используется трекером как сигнал к останову интегрирования и записи траектории на диск.'),

  h2('5.4', 'Математическое описание'),
  body('Финишная плоскость задаётся тремя точками a, b, c в глобальной системе координат, не лежащими на одной прямой. Указанные точки однозначно определяют саму плоскость; на их основе алгоритм строит единичную нормаль и применяет её для ориентированного знакового теста.'),
  body('Единичная нормаль к плоскости определяется через векторное произведение двух её рёбер:'),
  mathPara('n̂ = ', frac(['(b − a) × (c − a)'], ['|(b − a) × (c − a)|'])),
  body('где:'),
  mdash('a,  b,  c', 'три задающие точки из таблицы experiments.csv, м.'),
  mdash('n̂', 'единичная нормаль к плоскости (модуль 1). Векторное произведение даёт вектор, ортогональный плоскости; деление на его длину обеспечивает единичную норму.'),
  body('Поскольку ориентация n̂, полученная из векторного произведения, произвольна, она переворачивается таким образом, чтобы частица начинала движение с «обратной» стороны плоскости. Условие переворота записывается в виде:'),
  mathPara('n̂ ← −n̂,    если (', sub('x', '0'), ' − a) · n̂ > 0'),
  body('где:'),
  mdash(sub('x', '0'), 'начальная позиция частицы, м (particle_groups.csv).'),
  body('Скалярное произведение в левой части представляет собой знаковую проекцию вектора x₀ − a на n̂. Таким образом, после данной коррекции частица по построению стартует со знакового расстояния меньше нуля.'),
  body('Знаковое расстояние произвольной точки x от плоскости определяется выражением:'),
  mathPara('d(x) = (x − a) · n̂'),
  body('где:'),
  mdash('d > 0', 'точка x находится на той же стороне плоскости, куда направлена n̂.'),
  mdash('d < 0', 'точка x расположена на противоположной стороне.'),
  mdash('|d|', 'равно евклидову расстоянию от x до ближайшей точки плоскости.'),
  body('Между двумя последовательными отсчётами x_{n−1} и x_n (со знаковыми расстояниями d_{n−1} и d_n соответственно) пересечение фиксируется при смене знака:'),
  mathPara('( ', sub('d', 'n−1'), ' < 0   ∧   ', sub('d', 'n'), ' ≥ 0 )    ∨    ( ',
           sub('d', 'n−1'), ' ≥ 0   ∧   ', sub('d', 'n'), ' < 0 )'),
  body('Следует отметить, что использование полуоткрытой формы (d_n ≥ 0 / d_n < 0) гарантирует, что вырожденный случай d_n = 0 не будет пропущен.'),

  h2('5.5', 'Алгоритм решения'),
  body('Инициализация (один раз на частицу):'),
  code(
    '1. u ← b − a;  v ← c − a;  n ← u × v             ; вектор, ортогональный плоскости\n' +
    '2. нормализовать n; при вырождении положить n ← (1, 0, 0)   ; единичная нормаль\n' +
    '3. если (startPos − a) · n > 0:  n ← −n          ; стартовая точка должна иметь d < 0\n' +
    '4. callCount ← 0;  lastD ← 0;  initialized ← false   ; стартовое состояние per-particle'),
  body('Запрос (на каждом шаге):'),
  code(
    '1. callCount ← callCount + 1                     ; счётчик обращений к детектору\n' +
    '2. если callCount ≤ skipSteps:  вернуть false   ; пропуск первых N вызовов\n' +
    '                                                ; (нужно при пересечении плоскости несколько раз)\n' +
    '3. d ← (pos − a) · n                            ; знаковое расстояние до плоскости\n' +
    '4. если not initialized:                        ; первый «нормальный» вызов после skipSteps:\n' +
    '       lastD ← d; initialized ← true; вернуть false   ; запомнить d, без проверки\n' +
    '5. crossed ← (lastD < 0 ∧ d ≥ 0) ∨ (lastD ≥ 0 ∧ d < 0)   ; смена знака → пересечение\n' +
    '6. lastD ← d                                    ; сохранить для следующего шага\n' +
    '7. вернуть crossed'),
];

// ---------------------------------------------------------------------------
// 6. Distance to polyline
// ---------------------------------------------------------------------------
startSection('6');
const alg6 = [
  h1('6', 'Расстояние от точки до полилинии проектной орбиты'),
  body('Исходный код: src/optimizer.cpp:66, функция distanceToPolyline.'),

  h2('6.1', 'Назначение и характеристика'),
  body('Функция distanceToPolyline применяется для вычисления минимального евклидова расстояния от точки до ломаной, заданной упорядоченным списком вершин. Данная величина используется в качестве амплитуды огибающей в envelope-компоненте функции потерь (см. раздел 7), а также в скрипте analyze_trajectory_stability.py — для оценки отклонения траектории частицы от проектной орбиты.'),

  h2('6.2', 'Используемая информация'),
  tableCaption('6.1', 'Входные данные алгоритма'),
  table2col(
    ['Обозначение', 'Смысл'],
    [
      ['p', 'точка опроса в глобальных координатах'],
      ['poly', 'упорядоченный список вершин полилинии (N ≥ 2)'],
    ]),

  h2('6.3', 'Результаты решения'),
  body('Результатом работы функции является скалярное расстояние, удовлетворяющее условию ≥ 0. В случае, когда число вершин полилинии меньше двух, функция возвращает 0.'),

  h2('6.4', 'Математическое описание'),
  body('Полилиния представляет собой последовательность отрезков. Расстояние от точки до полилинии определяется как минимум расстояний до отдельных отрезков; задача проекции точки на отрезок решается стандартным методом — через клиппинг параметра проекции.'),
  body('Для отрезка [a, b] точка на бесконечной прямой, ближайшая к p, имеет параметр, определяемый выражением:'),
  mathPara(sup('t', '*'), ' = ', frac('(p − a) · (b − a)', sup('|b − a|', '2'))),
  body('где:'),
  mdash('a,  b', 'концы отрезка, м.'),
  mdash('p', 'точка опроса, м.'),
  mdash(sup('t', '*'), 'параметр прямой; значение t* = 0 соответствует a, t* = 1 — b; для t* ∈ (0, 1) проекция попадает на сам отрезок, а при выходе за [0, 1] — на его продолжение.'),
  mdash(sup('|b − a|', '2'), 'квадрат длины отрезка. Случай |b − a|² ≤ 10⁻³⁰ обрабатывается отдельно — это позволяет корректно обработать вырожденные отрезки нулевой длины.'),
  body('Чтобы ограничить проекцию самим отрезком (а не его продолжением), параметр t* клиппится в интервал [0, 1]:'),
  mathPara('t = clamp(', sup('t', '*'), ', 0, 1) = min(1,   max(0, ', sup('t', '*'), '))'),
  body('Ближайшая к p точка на отрезке и расстояние до неё вычисляются по формулам:'),
  mathPara('c = a + t (b − a)'),
  mathPara(sub('d', 'i'), ' = |p − c|'),
  body('где:'),
  mdash('c', 'основание перпендикулярной проекции; в случае выхода перпендикуляра за отрезок — ближайший к p его конец.'),
  mdash(sub('d', 'i'), 'евклидово расстояние от p до отрезка с индексом i, м.'),
  body('Таким образом, расстояние до полилинии определяется как минимум по всем N − 1 отрезкам:'),
  mathPara('dist(p, Π) = ', sub('min', 'i = 0..N−2'), '   ', sub('d', 'i'),
           ' = ', rad([sub('min', 'i'), '  ', sup(['|p − ', sub('c', 'i'), '|'], '2')])),
  body('где:'),
  mdash(['Π = (', sub('p', '0'), ', ', sub('p', '1'), ', …, ', sub('p', 'N−1'), ')'], 'полилиния, представленная упорядоченным списком из N вершин.'),
  body('Стоит отметить, что извлечение квадратного корня производится однократно — на финальном шаге. Внутри основного цикла сравниваются квадраты расстояний, что одновременно ускоряет вычисления и устраняет лишнюю операцию корня в горячем цикле.'),

  h2('6.5', 'Алгоритм решения'),
  code(
    'Вход:  p, poly[0..N-1]\n' +
    'Выход: расстояние dist ≥ 0\n\n' +
    '1. если N < 2:  вернуть 0                      ; полилиния не определена\n' +
    '2. best² ← +∞                                  ; текущий минимум квадрата расстояния\n' +
    '3. для i из [0..N-2]:                          ; перебор всех N−1 отрезков\n' +
    '     a ← poly[i];  b ← poly[i+1]               ; концы текущего отрезка\n' +
    '     u ← b − a;  len² ← u·u                    ; направляющий вектор и его |·|²\n' +
    '     t ← clamp( (p − a)·u / len²,  0,  1 )     ; параметр проекции с клиппингом в [0,1]\n' +
    '         если len² > 1e-30,  иначе 0           ; защита от вырожденного отрезка\n' +
    '     c ← a + t·u                               ; ближайшая к p точка отрезка\n' +
    '     d² ← (p − c)·(p − c)                      ; квадрат расстояния (без sqrt в цикле)\n' +
    '     если d² < best²: best² ← d²               ; обновление минимума\n' +
    '4. вернуть √best²                              ; sqrt извлекается однократно — в конце'),
];

// ---------------------------------------------------------------------------
// 7. Composite loss
// ---------------------------------------------------------------------------
startSection('7');
const alg7 = [
  h1('7', 'Составная функция потерь оптимизатора'),
  body('Исходный код: src/optimizer.cpp, функция runTrackingAndLoss.'),

  h2('7.1', 'Назначение и характеристика'),
  body('Составная функция потерь представляет собой скалярную целевую функцию, минимизируемую оптимизатором токов магнитов. Функция складывается из трёх аддитивных компонентов, отвечающих за различные физические дефекты пучка:'),
  dash('position_loss — обеспечивает прохождение частиц через центры магнитов;'),
  dash('direction_loss — обеспечивает параллельность выходящего пучка заданному направлению;'),
  dash('envelope_loss — препятствует росту бетатронных амплитуд.'),
  body('Стоит отметить, что компонент position_loss обновляется инкрементально в цикле трекинга — без сохранения траектории. Direction- и envelope-компоненты, напротив, вычисляются однократно на финальном шаге каждой частицы.'),

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
  body('Полная целевая функция представляет собой сумму трёх слагаемых:'),
  mathPara('L = ', sub('L', 'pos'), ' + ', sub('L', 'dir'), ' + ', sub('L', 'env')),
  body('где:'),
  mdash('L', 'скалярная потеря, минимизируемая оптимизатором. Единицы измерения — м² до взвешивания; коэффициенты веса нормируют слагаемые до сопоставимых порядков величин.'),
  body('Компонент Position loss отвечает за прохождение частиц через центры магнитов (а также заданных пользователем виртуальных целей) и вычисляется по формуле:'),
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
  body('где:'),
  mdash(sub('N', 'p'), 'число трекируемых частиц (5: опорная частица и 4 спутника).'),
  mdash('M', 'множество учитываемых в потере магнитов и виртуальных целей, формируемое после фильтрации по параметрам magnet_sp_numbers и exclude_device_types из JSON-конфига оптимизатора.'),
  mdash(sub('c', 'm'), 'геометрический центр магнита m в глобальных координатах, м.'),
  mdash(sub('w', 'm'), 'индивидуальный вес цели (по умолчанию принимается равным 1). Виртуальные цели из extra_target_magnets имеют, как правило, бо́льшие веса — порядка 20 или 50, что обеспечивает закрепление финальной сходимости.'),
  mdash([sub('x', 'i'), '(t)'], 'позиция частицы i на шаге интегрирования t, м.'),
  mdash([sub('min', 't'), '  ', sup(['|', sub('x', 'i'), '(t) − ', sub('c', 'm'), '|'], '2')], 'минимальный квадрат расстояния от траектории частицы i до центра магнита m, взятый по всей траектории. Величина вычисляется инкрементально: хранится только текущий минимум, а сама траектория не сохраняется.'),
  body('Компонент Direction loss обеспечивает параллельность импульсов на финишной плоскости заданному эталонному направлению и определяется выражением:'),
  mathPara(sub('L', 'dir'), ' = ', sub('w', 'dir'), '  ',
           new MathSum({
             children: [sup('sin', '2'), '(', sub('θ', 'i'), ')'].map(_m),
             subScript: _arr('i = 1..N_p'),
           })),
  mathPara('cos ', sub('θ', 'i'), ' = ',
           frac([subSup('p', 'i', 'final'), ' · r̂'], ['|', subSup('p', 'i', 'final'), '|'])),
  body('где:'),
  mdash(sub('w', 'dir'), 'весовой коэффициент направленной компоненты. Типичное значение — порядка 10³; столь большая величина обусловлена тем, что sin²(θ) безразмерен и обычно принимает малые значения.'),
  mdash(subSup('p', 'i', 'final'), 'импульс частицы i в момент пересечения финишной плоскости, кг·м/с.'),
  mdash('r̂', 'единичный эталонный вектор из конфига оптимизатора (direction_loss.reference). Как правило, в качестве эталона принимается касательная к проектной орбите в точке финиша.'),
  mdash(sub('θ', 'i'), 'угол между импульсом p_i^final и эталонным вектором r̂, рад.'),
  mdash([sup('sin', '2'), ' ', sub('θ', 'i'), ' = 1 − ', sup('cos', '2'), ' ', sub('θ', 'i')], 'данная величина квадратично зависит от θ при малых углах, благодаря чему компонент потерь принимает форму least-squares угловой ошибки.'),
  body('Компонент Envelope loss препятствует росту бетатронных колебаний при пролёте частиц через решётку и записывается в виде:'),
  mathPara(sub('L', 'env'), ' = ', sub('w', 'env'), '  ',
           new MathSum({
             children: [sup(['max(0,  ', subSup('A', 'i', 'final'), ' − ', subSup('A', 'i', 'init'), ')'], '2')].map(_m),
             subScript: _arr('i = 1..N_p'),
           })),
  mathPara(sub('A', 'i'), '(t) = dist( ', sub('x', 'i'), '(t), ', sub('Π', 'design'), ' )'),
  body('где:'),
  mdash(sub('w', 'env'), 'весовой коэффициент огибающей. Типичное значение составляет порядка 10⁶; столь большая величина связана с тем, что амплитуда A измеряется в метрах, и даже её метровый рост является недопустимым.'),
  mdash([sub('A', 'i'), '(t)'], 'амплитуда огибающей частицы i — её расстояние от проектной полилинии в момент времени t. Величина рассчитывается по алгоритму, описанному в § 6.'),
  mdash(sub('Π', 'design'), 'проектная полилиния (центры элементов решётки плюс прелюдийный/постлюдийный сегменты по касательной к орбите). Перестраивается в начале каждого прогона оптимизации.'),
  mdash([subSup('A', 'i', 'init'), ',  ', subSup('A', 'i', 'final')], 'амплитуды огибающей в начале и в конце прогона, м.'),
  body('Стоит подчеркнуть, что обёртка max(0, ·) означает следующее: «затухание» — то есть приближение пучка к орбите (когда A_i^final < A_i^init) — не штрафуется. Штраф накладывается только на рост амплитуды.'),

  h2('7.5', 'Алгоритм решения'),
  body('А. Подготовка статических данных (выполняется один раз за вызов):'),
  code(
    '1. magnets_idx ← фильтрация элементов по             ; учитываемые магниты: только нужные СП,\n' +
    '       (sp_number ∈ S) ∧ (device_type ∉ exclude)     ; кроме исключённых типов\n' +
    '2. virtual_targets ← список (c, w) из extra_target_magnets   ; виртуальные цели с весами\n' +
    '3. design_poly ← центры элементов 1..N +             ; проектная орбита:\n' +
    '       прелюдия  (экстраполяция первого сегмента на −2.94 м) +\n' +
    '       постлюдия (экстраполяция последнего сегмента вперёд)'),
  body('Б. Цикл по каждой частице i:'),
  code(
    'minDist²[m] ← +∞ для каждого m из magnets_idx ∪ virtual_targets   ; инкрементальный аккум.\n' +
    'A_init ← distanceToPolyline(x_init, design_poly)    ; начальная амплитуда огибающей\n' +
    'для каждого шага до пересечения финишной плоскости:\n' +
    '    шаг Бориса (borisStep)                          ; продвижение частицы (см. § 1)\n' +
    '    для каждого m из бакета:                        ; учёт только близких магнитов\n' +
    '        d² ← |x − c_m|²                             ; квадрат расстояния до центра магнита\n' +
    '        если d² < minDist²[m]: minDist²[m] ← d²    ; обновление инкрементального минимума\n' +
    'A_final ← distanceToPolyline(x_final, design_poly)  ; финальная амплитуда огибающей\n' +
    'θ_i     ← угол(p_final, r̂)                          ; угол между p_final и эталоном'),
  body('В. Агрегация:'),
  code(
    'position_loss  = Σ_(i,m) w_m · minDist²[i, m]            ; § 7.4: близость к магнитам\n' +
    'direction_loss = w_dir · Σ_i sin²(θ_i)                   ; параллельность пучка эталону\n' +
    'envelope_loss  = w_env · Σ_i max(0, A_final − A_init)²   ; штраф только за РОСТ огибающей\n' +
    'loss           = position_loss + direction_loss + envelope_loss'),
  body('Если хотя бы одна частица не достигает финишной плоскости за 50 миллионов шагов, функция возвращает признак сбоя, и оптимизатор отвергает текущее состояние.'),
];

// ---------------------------------------------------------------------------
// 8. Stochastic optimizer
// ---------------------------------------------------------------------------
startSection('8');
const alg8 = [
  h1('8', 'Стохастический оптимизатор (hill-climber c киками)'),
  body('Исходный код: src/optimizer.cpp, функция main и вспомогательные функции.'),

  h2('8.1', 'Назначение и характеристика'),
  body('Оптимизатор представляет собой однопоточный жадный hill-climber с периодическими «киками» (kicks) по всему вектору параметров, применяемыми для выхода из локальных минимумов. В качестве управляемых параметров рассматриваются токи магнитов (arg_val) и индивидуальные калибровочные множители (correct_coeff), перечисленные в секции parameter_groups_per_position JSON-конфига.'),
  body('Реализация обладает следующими свойствами:'),
  dash('асимметричный шаг — на каждой итерации одна случайная группа возмущается на ±step_abs (для тока) либо ±step_percent · current (для калибровочного множителя);'),
  dash('кики — каждые kick_every итераций все группы получают возмущение ±kick_scale одновременно;'),
  dash('жадный приём — новое состояние принимается только при уменьшении функции потерь; в противном случае выполняется откат;'),
  dash('контрольные точки — полная история сохраняется каждые CHECKPOINT_EVERY итераций, что обеспечивает возможность восстановления после прерывания.'),

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
  body('Оптимизатор выполняет обход пространства параметров методом Монте-Карло с жадным приёмом. Прежде всего введём следующие обозначения:'),
  mdash(['θ ∈ ', sup('ℝ', 'M')], 'вектор всех изменяемых параметров (токи магнитов и индивидуальные коэффициенты), имеющий длину M.'),
  mdash('G', 'число групп, объявленных в parameter_groups_per_position. Каждая группа g может объединять несколько слотов параметров; мутация группы g изменяет все её слоты одновременно.'),
  mdash('L(θ)', 'составная функция потерь, определённая в § 7.'),
  mdash([sub('θ', 'best')], 'текущий лучший вектор параметров; соответствующее значение функции потерь обозначается как L_best = L(θ_best).'),
  body('Обычная итерация состоит в следующем: случайным образом выбирается одна группа g, и её параметры сдвигаются на индивидуальный шаг со случайным знаком. Преобразование записывается в виде:'),
  mathPara(sup(sub('θ', 'g'), '′'), ' = ', sub('θ', 'g'), ' + ', sub('ε', 'g'), ' · sign(u − ½)'),
  mathPara(sup(sub('θ', 'k'), '′'), ' = ', sub('θ', 'k'), ',    для всех k ≠ g'),
  body('где:'),
  mdash('u ~ Uniform(0, 1)', 'случайное число, генерируемое заново на каждой итерации.'),
  mdash('sign(u − ½) ∈ {−1, +1}', 'знак предложенного сдвига; положительное и отрицательное направления равновероятны.'),
  mdash(sub('ε', 'g'), 'размер шага для группы g, определяемый типом колонки изменяемого параметра:'),
  mathPara(sub('ε', 'g'), ' = ', sub('step_abs', 'g'), ',    если группа g целит колонку arg_val'),
  mathPara(sub('ε', 'g'), ' = ', sub('step_percent', 'g'), '  · ', sub('θ', 'g'), ',    если группа g целит колонку correct_coeff'),
  body('где:'),
  mdash(sub('step_abs', 'g'), 'абсолютный шаг в амперах; типичные значения лежат в диапазоне 0,1–1,0.'),
  mdash(sub('step_percent', 'g'), 'относительный шаг (типично 0,05, что соответствует 5%), пропорциональный текущему значению. Такой выбор обеспечивает осмысленное относительное возмущение даже при уменьшении θ_g.'),
  body('Кик-итерация выполняется каждые kick_every итераций. На ней все группы возмущаются одновременно, что позволяет покидать локальные минимумы:'),
  mathPara(sup(sub('θ', 'g'), '′'), ' = ', sub('θ', 'g'), ' + ', sub('δ', 'g'), ',    для всех g'),
  mathPara(sub('δ', 'g'), ' ~ Uniform(−Δ, +Δ)'),
  body('где:'),
  mdash('Δ = kick_scale', 'амплитуда кика (амперы или доля от текущего значения — в зависимости от типа колонки группы).'),
  body('Стоит отметить, что независимые равномерные выборки означают независимое движение каждой группы. При этом совместное предложение оказывается достаточно крупным, чтобы покинуть тот бассейн притяжения, который пошаговый поиск никогда бы не покинул.'),
  body('Жадный приём, общий для обоих типов итераций, формулируется следующим образом: предложение принимается только при условии уменьшения лучшего значения функции потерь:'),
  mathPara('θ,  ', sub('L', 'best'), '  ←  θ′,  L(θ′),    если L(θ′) < ', sub('L', 'best')),
  mathPara('θ,  ', sub('L', 'best'), '  ←  θ,  ', sub('L', 'best'), ',    иначе (откат)'),
  body('Таким образом, рассматриваемый алгоритм представляет собой чистый hill-climber — без введения «температуры» и без Metropolis-приёма ухудшающих состояний. Возможность покинуть локальный минимум обеспечивается исключительно за счёт кик-итераций.'),

  h2('8.5', 'Алгоритм решения'),
  code(
    '1. Загрузить конфигурацию:                       ; стадия инициализации (~1 с)\n' +
    '     - прочитать JSON, заполнить OptimConfig\n' +
    '     - прочитать sp_elements, hardware_library, lattice_configs,\n' +
    '       particle_groups                            ; вся БД — в память один раз\n' +
    '     - отфильтровать строки решётки по lattice_id и lattice_max_gid\n' +
    '                                                  ; для экспериментов на части кольца\n' +
    '2. Построить ускоритель:\n' +
    '     - загрузить карты поля в общий кеш (один .bin на device_type)   ; mapCache\n' +
    '     - для каждой строки: создать MagneticElement, computeBounds()   ; bounding sphere\n' +
    '     - acc.buildSpatialIndex()                    ; bucket-индекс (см. § 3)\n' +
    '3. Построить список параметров:                   ; θ — изменяемые слоты\n' +
    '     - для каждой записи parameter_groups_per_position:\n' +
    '         найти соответствующие (sp_type, sp_eid)  ; какие элементы группа охватывает\n' +
    '         собрать индексы элементов\n' +
    '         сохранить ссылки на изменяемые поля     ; arg_val / correct_coeff in-place\n' +
    '4. Начальная оценка:\n' +
    '     L_best ← runTrackingAndLoss(θ)               ; baseline для жадного приёма\n' +
    '     при сбое — прервать                          ; стартовая конфигурация нерабочая\n' +
    '5. Цикл i = 1..n_iterations:                      ; ГОРЯЧИЙ ЦИКЛ оптимизации\n' +
    '     если i % kick_every == 0:                    ; КИК-ИТЕРАЦИЯ\n' +
    '         для каждой группы g:                     ; все группы одновременно\n' +
    '             θ_g += Uniform(−kick_scale, +kick_scale)\n' +
    '     иначе:                                       ; ОБЫЧНАЯ ИТЕРАЦИЯ\n' +
    '         выбрать случайную группу g\n' +
    '         θ_g_old ← θ_g                            ; сохранить для возможного отката\n' +
    '         θ_g ← θ_g + step · sign(rand − 0.5)     ; шаг со случайным знаком\n' +
    '     L_new ← runTrackingAndLoss(θ)                ; полный прогон трекинга\n' +
    '     если L_new — сбой или L_new ≥ L_best:        ; ЖАДНЫЙ ПРИЁМ\n' +
    '         откатить мутации этой итерации          ; неудачное предложение — назад\n' +
    '     иначе:\n' +
    '         L_best ← L_new                           ; принимаем улучшение\n' +
    '         добавить (i, L_best, snapshot θ) в историю   ; для viz_optim_history.py\n' +
    '     если i % CHECKPOINT_EVERY == 0:              ; периодический чекпойнт\n' +
    '         сохранить history.json (атомарно через временный файл)\n' +
    '6. Финальная фаза:                                ; результаты — на диск\n' +
    '     - записать лучшее θ обратно в sp_elements.csv и\n' +
    '       hardware_library.csv                       ; пользователь видит итог\n' +
    '     - пересобрать lattice_configs.csv (сохранить суб-решётки)\n' +
    '     - выполнить финальный прогон, записать trajectories_exp{id}.bin\n' +
    '                                                  ; чтобы visualize.py показал итог\n' +
    '     - сохранить итоговый history.json            ; с пометкой "final": true'),
];

// ---------------------------------------------------------------------------
// 9. LTE -> lattice
// ---------------------------------------------------------------------------
startSection('9');
const alg9 = [
  h1('9', 'Генерация решётки из LTE-описания'),
  body('Исходный код: lte_to_lattice.py.'),

  h2('9.1', 'Назначение и характеристика'),
  body('Скрипт lte_to_lattice.py выполняет разбор файлов формата Elegant LTE (SYLA_desc.lte, SYLA_line.lte), классифицирует элементы согласно конфигурации lte_mapping.json и производит развёртку кольца вдоль хорды. Результатом работы являются редактируемые шаблоны sp_elements.csv и геометрическая развёртка lattice_layout.csv.'),
  body('Реализация обладает следующими особенностями:'),
  dash('секционирование DL/JL — длинные изгибающие магниты разбиваются на 5 фрагментов, размещаемых на хорде вместо дуги;'),
  dash('сохранение пользовательских значений arg_val — при повторных запусках значения, ранее введённые в sp_elements.csv, не затираются.'),

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
  body('Длина хорды поворотного магнита вычисляется следующим образом. Для магнита с заданной длиной дуги L и углом поворота α радиус кривизны определяется как ρ = L/α. Длина прямой хорды, соединяющей входную и выходную грани дуги, задаётся выражением:'),
  mathPara(sub('L', 'chord'), ' = 2 ρ sin(α/2) = ',
           frac('2 L sin(α/2)', 'α')),
  body('где:'),
  mdash('L', 'длина дуги поворотного магнита, м; берётся из параметра L= LTE-описания.'),
  mdash('α', 'полный угол поворота, рад; берётся из параметра ANGLE= LTE. Положительное значение α соответствует повороту влево, отрицательное — вправо.'),
  mdash('ρ = L/α', 'радиус кривизны магнита, м.'),
  mdash(sub('L', 'chord'), 'длина хорды, по которой выпрямляется магнит. При α ≠ 0 хорда оказывается чуть короче дуги L; при α → 0 хорда стремится к L, как и ожидается.'),
  body('Магнит заменяется пятью тонкими ломтиками, равномерно расположенными по хорде в параметрических точках τ ∈ {0, 0,25, 0,5, 0,75, 1,0}. Каждый ломтик ориентирован по направлению хорды.'),
  body('Развёртка кольца производится сегмент за сегментом. Каждый LTE-элемент длины ℓ задаёт в плоскости поворот Δψ (значение Δψ = 0 соответствует дрейфу, Δψ = α — поворотному магниту) и сдвиг геометрической середины. Обновление состояния записывается в виде:'),
  mathPara('ψ′ = ψ + Δψ'),
  mathPara('x′ = x + ℓ cos( ψ + ', frac('Δψ', '2'), ' )'),
  mathPara('y′ = y + ℓ sin( ψ + ', frac('Δψ', '2'), ' )'),
  body('где:'),
  mdash('ψ', 'yaw-угол — поворот вокруг глобальной оси z, рад.'),
  mdash('Δψ', 'yaw-поворот, наработанный текущим элементом, рад. Для поворотного магнита Δψ = α, для прямого элемента Δψ = 0.'),
  mdash('ℓ', 'длина элемента по хорде; для поворотного магнита соответствует пятой части L_chord, приходящейся на один ломтик.'),
  mdash('x,  y', 'накопленные координаты текущей позиции в глобальной системе, м.'),
  body('Стоит отметить, что использование «среднего» значения yaw ψ + Δψ/2 в тригонометрических аргументах — вместо входного ψ или выходного ψ + Δψ — позволяет сохранить касательность полилинии истинной дуге на концах ломтика и симметрично уравновесить ошибку кривизны вдоль него.'),

  h2('9.5', 'Алгоритм решения'),
  code(
    '1. Чтение JSON-конфига → mapping, regex-шаблоны, пороги\n' +
    '                                                ; правила соответствия LTE-типов и СП\n' +
    '2. Разбор SYLA_desc.lte:                        ; БИБЛИОТЕКА УСТРОЙСТВ\n' +
    '     для каждой строки "ИМЯ: ТИП, key1=val1, ..."\n' +
    '         классифицировать ТИП через конфиг      ; DL/JL/SBM/DP/DQ/QP/SQ/COR/MARK/...\n' +
    '         сохранить (ИМЯ, классиф_тип, длина, угол, K1, ...)\n' +
    '3. Разбор SYLA_line.lte:                        ; ПОСЛЕДОВАТЕЛЬНОСТЬ ЭЛЕМЕНТОВ\n' +
    '     развернуть вложенные LINE-выражения        ; LTE поддерживает рекурсивные ссылки\n' +
    '     в плоскую последовательность (ИМЯ, повтор)\n' +
    '4. Классификация в (sp_type, sp_element_id):    ; разбиение на типовые СП\n' +
    '     standard — обычные СП, повторяющиеся 1..N\n' +
    '     junction — стыковочные СП 18-21 (отличающаяся геометрия)\n' +
    '5. Развёртка кольца:                            ; РАЗМЕЩЕНИЕ В ПРОСТРАНСТВЕ\n' +
    '     для каждого элемента в плоской последовательности:\n' +
    '         если DL/JL: выдать 5 фрагментов на хорде   ; спрямление длинных бендов\n' +
    '         иначе:      выдать 1 фрагмент\n' +
    '         обновить (x, y, yaw) по формулам п. 9.4    ; продвижение курсора\n' +
    '6. Прочитать существующий sp_elements.csv (если есть):\n' +
    '     сохранить arg_val для совпадающих ключей    ; пользователь не теряет токи\n' +
    '     (sp_type, sp_eid)                           ; от прошлой оптимизации\n' +
    '7. Записать sp_elements.csv                     ; шаблоны, по одной строке на ключ\n' +
    '8. Записать lattice_layout.csv                  ; полная развёртка кольца'),
];

// ---------------------------------------------------------------------------
// 10. Sub-lattice assembly
// ---------------------------------------------------------------------------
startSection('10');
const alg10 = [
  h1('10', 'Сборка таблицы конфигурации и суб-решёток'),
  body('Исходный код: assemble_lattice_configs.py.'),

  h2('10.1', 'Назначение и характеристика'),
  body('Скрипт assemble_lattice_configs.py выполняет объединение редактируемого шаблона sp_elements.csv с геометрической развёрткой lattice_layout.csv в итоговую входную таблицу симулятора lattice_configs.csv. Дополнительно генерируются укороченные копии решётки с идентификаторами 102, 103, 104 и далее — на основании sub_lattices.json, — содержащие только первые max_gid строк основной решётки. Таким образом обеспечивается возможность запускать эксперименты на части кольца без повторной развёртки LTE-описаний.'),

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
  body('Стоит отметить, что какая-либо арифметика в данном алгоритме отсутствует. Вместо неё выполняется реляционное соединение (join) таблицы геометрии lattice_layout с таблицей токов и калибровок sp_elements, дополненное объединением с укороченными копиями основной решётки.'),
  body('В нотации построения множеств соответствующие соотношения записываются в виде:'),
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
  body('где:'),
  mdash([sub('L', '0'), ' = 101'], 'идентификатор основной решётки, задаваемый в lte_mapping.json.'),
  mdash('ℓ', 'отдельная строка таблицы lattice_layout.csv.'),
  mdash('ℓ.q = (ℓ.x, ℓ.y, ℓ.z, ℓ.yaw, ℓ.pitch, ℓ.roll)', 'геометрический кортеж элемента.'),
  mdash('s', 'в первом блоке формул — строка sp_elements.csv; во втором — дескриптор суб-решётки из sub_lattices.json.'),
  body('Условие соединения формулируется по составному ключу (sp_type, sp_element_id): каждая строка layout получает свой device_type_id, ток и калибровочный коэффициент из соответствующего SP-шаблона. Объединение по суб-решёткам создаёт дополнительные строки с новым идентификатором lattice_id, однако только для instance_id ≤ max_gid — то есть для префикса основного кольца.'),

  h2('10.5', 'Алгоритм решения'),
  code(
    '1. sp_dict ← load_sp_elements()              ; шаблоны СП по ключу (sp_type, sp_eid)\n' +
    '                                              ; → (device_type, arg_val, correct_coeff)\n' +
    '2. output_rows ← []                          ; накопитель результата\n' +
    '3. для каждой строки r из lattice_layout.csv:    ; РЕЛЯЦИОННЫЙ JOIN: геометрия + шаблон\n' +
    '     key ← (r.sp_type, r.sp_element_id)      ; составной ключ соединения\n' +
    '     если key не в sp_dict: выдать warning и пропустить\n' +
    '                                              ; элемент layout без шаблона — диагностика\n' +
    '     dev, arg, coeff ← sp_dict[key]          ; берём данные из SP-шаблона\n' +
    '     добавить запись {LATTICE_ID, r.global_id, dev, r.x..roll, arg, coeff}\n' +
    '                                              ; полная итоговая строка для симулятора\n' +
    '4. для каждого (sub_id, max_gid) из SUB_LATTICES:    ; ПОДРЕШЁТКИ для частичных экспериментов\n' +
    '     для r из output_rows, где r.instance_id ≤ max_gid:\n' +
    '         добавить копию с lattice_id ← sub_id    ; префикс кольца с новым идентификатором\n' +
    '5. Запись CSV с повторными попытками при PermissionError\n' +
    '   (блокировка файла в Windows — например, Excel держит файл открытым)'),
];

// ---------------------------------------------------------------------------
// 11. Envelope + FFT
// ---------------------------------------------------------------------------
startSection('11');
const alg11 = [
  h1('11', 'Анализ огибающей и бетатронных частот по FFT'),
  body('Исходный код: analyze_trajectory_stability.py.'),

  h2('11.1', 'Назначение и характеристика'),
  body('Скрипт analyze_trajectory_stability.py выполняет постобработку бинарного файла траектории для оценки динамической устойчивости пучка. В рамках анализа решаются следующие задачи:'),
  dash('построение огибающей A_i(t) = dist(x_i(t), design_polyline) для каждой из пяти частиц;'),
  dash('подгонка линейной и экспоненциальной моделей роста A_i(t) с оценкой коэффициента детерминации R²;'),
  dash('FFT-анализ A_i(t) для определения частоты бетатронных колебаний (tune);'),
  dash('формирование вердикта по каждой частице: STABLE, OSCILLATING, GROWING или DRIFTED.'),

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
  body('Амплитуда огибающей для каждой частицы i и каждого сохранённого отсчёта траектории k определяется выражением:'),
  mathPara(sub('A', 'i'), '(', sub('t', 'k'), ') = dist( ', sub('x', 'i'), '(', sub('t', 'k'), '),  ', sub('Π', 'design'), ' )'),
  body('где:'),
  mdash([sub('x', 'i'), '(', sub('t', 'k'), ')'], 'позиция частицы i в k-й сохранённой временной точке, м.'),
  mdash(sub('Π', 'design'), 'проектная полилиния — центры элементов решётки, дополненные прелюдийным и постлюдийным сегментами по касательной.'),
  mdash([sub('A', 'i'), '(', sub('t', 'k'), ')'], 'мгновенное евклидово расстояние от траектории до проектной орбиты, м; вычисляется по алгоритму § 6.'),
  body('Линейная подгонка реализует базовую дрейф-модель и записывается в виде:'),
  mathPara('A(t) ≈ m t + ', sub('A', '0')),
  mathPara('(m, ', sub('A', '0'), ') = argmin   ',
           new MathSum({
             children: [sup(paren(['A(', sub('t', 'k'), ') − m  ', sub('t', 'k'), ' − ', sub('A', '0')]), '2')].map(_m),
             subScript: _arr('k'),
           })),
  body('где:'),
  mdash('m', 'линейный темп роста огибающей, м/с. Положительное значение m свидетельствует о дрейфе пучка от проектной орбиты.'),
  mdash(sub('A', '0'), 'экстраполированное значение огибающей при t = 0, м.'),
  body('Стоит отметить, что минимизация выполняется стандартным методом наименьших квадратов — в замкнутой форме через вызов numpy.polyfit(..., deg=1).'),
  body('Экспоненциальная подгонка применяется для описания модели роста или затухания амплитуды и задаётся выражением:'),
  mathPara('A(t) ≈ ', sub('A', '0'), ' · ', sup('e', 'γ t')),
  mathPara('ln A(t) ≈ ln ', sub('A', '0'), ' + γ t'),
  body('где:'),
  mdash('γ', 'экспоненциальный темп роста огибающей, с⁻¹ — своего рода Ляпунов-подобный показатель. При γ > 0 амплитуда растёт, при γ < 0 — затухает.'),
  body('Следует отметить, что подгонка линейна по величине ln A. Соответственно, она выполняется в логарифмическом пространстве методом наименьших квадратов по точкам, удовлетворяющим условию A > A_floor (что позволяет избежать логарифма нуля).'),
  body('Качество экспоненциальной подгонки оценивается коэффициентом детерминации, определяемым выражением:'),
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
  body('где:'),
  mdash(sub('SS', 'res'), 'сумма квадратов остатков, характеризующая ошибку подгонки.'),
  mdash(sub('SS', 'tot'), 'суммарная сумма квадратов, характеризующая вариацию данных.'),
  mdash('⟨ln A⟩', 'среднее значение логарифма амплитуды по выборке.'),
  mdash([sup('R', '2'), ' ∈ (−∞, 1]'], 'коэффициент детерминации. Значения, близкие к единице, соответствуют лучшей подгонке; отрицательные значения свидетельствуют о том, что подгонка хуже, чем простое использование среднего.'),
  body('Таким образом, низкое значение R² дисквалифицирует оценку γ, даже если она получается численно большой.'),
  body('Бетатронная частота определяется через дискретное преобразование Фурье огибающей, выполняемое при равномерном шаге Δt_sample = t_k − t_{k−1}:'),
  mathPara('Â(', sub('f', 'n'), ') = ',
           new MathSum({
             children: ['A(', sub('t', 'k'), ') · ', sup('e', ['−2π i  ', sub('f', 'n'), '  ', sub('t', 'k')])].map(_m),
             subScript: _arr('k = 0..K−1'),
           })),
  mathPara(sub('f', 'n'), ' = ', frac('n', ['K  ', sub('Δt', 'sample')]), ',    n = 0, 1, …, K−1'),
  body('где:'),
  mdash('K', 'число отсчётов огибающей в траектории.'),
  mdash(sub('f', 'n'), 'дискретный частотный бин, Гц.'),
  mdash(['Â(', sub('f', 'n'), ')'], 'комплексный коэффициент Фурье на частоте f_n.'),
  body('Доминирующая бетатронная частота определяется как частота того бина, который имеет максимальную амплитуду, за исключением DC-бина (содержащего только среднее значение сигнала):'),
  mathPara(sub('f', 'betatron'), ' = ', sub('argmax', [sub('f', 'n'), ' > ', sub('f', 'min')]),
           '   |Â(', sub('f', 'n'), ')|'),
  body('где:'),
  mdash(sub('f', 'min'), 'малая пороговая частота, маскирующая DC-компоненту, а также (при необходимости) несколько соседних бинов — для устранения утечки от ненулевого среднего значения сигнала.'),
  mdash(sub('f', 'betatron'), 'частота наиболее выраженного спектрального пика, Гц. Для перевода полученной величины в дробный «тюн» (tune) её необходимо разделить на частоту обращения по кольцу; анализатор выполняет данное преобразование автоматически, когда длина решётки известна.'),

  h2('11.5', 'Алгоритм решения'),
  code(
    '1. Загрузить траектории (5 массивов точек (Ni, 3)).      ; trajectories_exp{id}.bin\n' +
    '2. Определить lattice_id по experiments.csv, загрузить lattice_configs,\n' +
    '   построить design_poly (центры + прелюдия/постлюдия).  ; проектная орбита\n' +
    '3. Для каждой частицы i:\n' +
    '     a. A_i ← [distanceToPolyline(p, design_poly) для p из traj_i]\n' +
    '                                                         ; огибающая — по алгоритму § 6\n' +
    '     b. t_i ← равномерный шаг (dt × traj_stride × stride) ; восстановление времени отсчётов\n' +
    '     c. линейная подгонка: m, A₀ ← polyfit(t_i, A_i, deg=1)   ; baseline дрейф-модель\n' +
    '     d. экспоненциальная подгонка (только точки с A > A_floor):   ; защита от log 0\n' +
    '            log A = γ·t + log A₀; вычислить R²                   ; γ — Ляпунов-показатель\n' +
    '     e. FFT(A_i); маскировать DC-бин;                     ; убираем влияние среднего\n' +
    '        выбрать максимум как f_betatron                   ; доминирующая частота колебаний\n' +
    '     f. классификация (вердикт по частице):\n' +
    '            DRIFTED      если A_init ≈ 0 и m_lin > δ     ; растущий дрейф от орбиты\n' +
    '            GROWING      если γ > γ_threshold и R² ≥ R²_min   ; экспоненциальный рост\n' +
    '            OSCILLATING  если есть выраженный пик FFT и γ ≤ γ_threshold\n' +
    '                                                         ; устойчивые бетатронные колебания\n' +
    '            STABLE       иначе                            ; никаких аномалий\n' +
    '4. Построить графики, вывести таблицу вердиктов.'),
];

// ---------------------------------------------------------------------------
// Приложения
// ---------------------------------------------------------------------------
startSection(null);
const appA = [
  h1Unnumbered('Приложение А. Файлы конфигурации', { id: 'A' }),
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
  h1Unnumbered('Приложение Б. Файлы данных', { id: 'B' }),
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
  ...TOC_ENTRIES.map(([anchor, num, text]) => tocItem(anchor, num, text)),
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
