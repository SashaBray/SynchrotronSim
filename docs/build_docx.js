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
  PageNumber, PageBreak, TabStopType, TabStopPosition
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
  body('Уравнения движения релятивистской заряженной частицы в электромагнитном поле:'),
  code('dx/dt = v\ndp/dt = q · (E + v × B)\np      = γ · m · v,    γ = √(1 + p²/(m²·c²))'),
  body('Шаг расщепляется на семь подопераций (drift-kick-drift):'),
  code(
    '1) Полудрейф:       x_{n+½} = x_n + (p_n / (γ_n · m)) · dt/2\n' +
    '2) Опрос поля:      B, E ← getTotalFields(x_{n+½})\n' +
    '3) Полу-E-толчок:   p⁻ = p_n + q · E · dt/2\n' +
    '4) γ⁻ = √(1 + |p⁻|² / (m²·c²))\n' +
    '5) Магнитный поворот:\n' +
    '       t  = q · B · dt / (2 · γ⁻ · m)\n' +
    '       s  = 2 · t / (1 + |t|²)\n' +
    '       p\' = p⁻ + (p⁻ × t)\n' +
    '       p⁺ = p⁻ + (p\' × s)            |p⁺| = |p⁻| строго\n' +
    '6) Полу-E-толчок:   p_{n+1} = p⁺ + q · E · dt/2\n' +
    '7) Полудрейф:       x_{n+1} = x_{n+½} + (p_{n+1} / (γ_{n+1} · m)) · dt/2'),

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
  body('Пусть (i, j, k) — нижние целые индексы узла, (tx, ty, tz) — дробные смещения в диапазоне [0, 1]. Трилинейная интерполяция:'),
  code(
    'F(p) =  v000·(1-tx)(1-ty)(1-tz) + v100·tx·(1-ty)(1-tz)\n' +
    '      + v010·(1-tx)·ty·(1-tz)   + v110·tx·ty·(1-tz)\n' +
    '      + v001·(1-tx)(1-ty)·tz    + v101·tx·(1-ty)·tz\n' +
    '      + v011·(1-tx)·ty·tz       + v111·tx·ty·tz'),
  body('Правило YZ-симметрии: если useSymmetry истинен и p.x < 0, опрос производится в точке (-p.x, p.y, p.z), после чего у результата инвертируется x-компонента: F.x ← -F.x (антисимметрия Bx; By, Bz симметричны).'),

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
  body('Проекция ограничивающих сфер магнитов на плоскость X-Y (кольцо лежит в z ≈ 0):'),
  code('bucketSize = max(1.0, 2 · max_i worldRadius_i)'),
  body('Такой размер бакета гарантирует, что круг проекции любого магнита укладывается в группу 2×2 бакетов, то есть регистрируется не более 4 раз. Координаты бакета по точке запроса:'),
  code(
    'bx = floor(x / bucketSize) - bxMin\n' +
    'by = floor(y / bucketSize) - byMin'),

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
  code(
    'если |globalPos − worldCenter|² > worldRadius²:  return       (прескрининг)\n' +
    'localPos = R⁻¹ · (globalPos − T)                  (мир → локал)\n' +
    'localF   = trilinear(map, localPos, useSymmetry)\n' +
    'если localF == 0:                                 return      (разреж. сетка)\n' +
    'globalF  = R · localF                             (локал → мир)\n' +
    'scale    = correctCoeff · fieldDirection · (currentScale если isActive)\n' +
    'totalB (или totalE) += scale · globalF             (по type)'),

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
  body('Нормаль к плоскости вычисляется по векторному произведению:'),
  code('n̂ = (b − a) × (c − a) / |(b − a) × (c − a)|'),
  body('Знак нормали выбирается так, чтобы начальная точка лежала на отрицательной стороне:'),
  code('n̂ ← −n̂,   если (startPos − a) · n̂ > 0'),
  body('Знаковое расстояние от произвольной точки pos:'),
  code('d(pos) = (pos − a) · n̂'),
  body('Условие пересечения между двумя последовательными отсчётами:'),
  code('(lastD < 0 ∧ d ≥ 0)  ∨  (lastD ≥ 0 ∧ d < 0)'),

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
  body('Для каждого сегмента [a, b] точка p проектируется на прямую, проходящую через сегмент. Параметр проекции t ограничивается интервалом [0, 1] (clamp):'),
  code(
    't = clamp( ((p − a) · (b − a)) / |b − a|²,  0,  1 )\n' +
    'c = a + t · (b − a)\n' +
    'd_i² = |p − c|²'),
  body('Итоговое расстояние:'),
  code('dist(p, poly) = √( min_i d_i² )'),

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
  code(
    'position_loss  = Σ_(i частица, m магнит) w_m · min_t |x_i(t) − c_m|²\n' +
    'direction_loss = w_dir · Σ_i sin²(θ_i),\n' +
    '                 cos θ_i = (p_i_final · r̂) / |p_i_final|\n' +
    'envelope_loss  = w_env · Σ_i max(0, A_i_final − A_i_init)²,\n' +
    '                 A_i(t) = distanceToPolyline(x_i(t), design_polyline)\n\n' +
    'loss = position_loss + direction_loss + envelope_loss'),
  body('Замечание: penality применяется только к росту огибающей; бетатронное затухание считается полезным эффектом и не штрафуется.'),

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
  body('Пусть θ ∈ ℝ^M — вектор управляемых параметров, L(θ) — функция потерь из раздела 7, g ∈ {1..G} — случайно выбранная группа.'),
  code(
    'Обычная итерация:\n' +
    '   θ\'_g = θ_g + ε_g · sign(rand() − 0.5),\n' +
    '   θ\'_k = θ_k для k ≠ g\n' +
    '   ε_g = step_abs          (если колонка группы = arg_val)\n' +
    '       или step_percent · θ_g (если колонка = correct_coeff)\n\n' +
    'Кик-итерация:\n' +
    '   θ\' = θ + δ,  δ_g ~ Uniform(−kick_scale, +kick_scale) ∀ g\n\n' +
    'Приём:\n' +
    '   θ ← θ\' если L(θ\') < L(θ_best); L(θ_best) ← L(θ\')\n' +
    '   иначе θ без изменения'),

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
  body('Для изгиба длины L с углом α длина хорды:'),
  code('Lc = 2 · (L/α) · sin(α/2)'),
  body('Пять фрагментов размещаются на долях 0, 0.25, 0.5, 0.75, 1.0 хорды с поворотом по yaw, соответствующим направлению хорды в данной точке.'),
  body('Перемещение по кольцу: для каждого сегмента длиной ℓ с поворотом Δψ выполняется обновление состояния:'),
  code(
    'yaw\' = yaw + Δψ\n' +
    'x\'   = x + ℓ · cos(yaw + Δψ/2)\n' +
    'y\'   = y + ℓ · sin(yaw + Δψ/2)'),

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
  body('Арифметика отсутствует — выполняется реляционное соединение:'),
  code(
    'output = {(LATTICE_ID, layout.global_id, sp.device_type_id,\n' +
    '           layout.x..roll, sp.arg_val, sp.correct_coeff)\n' +
    '         | layout ∈ lattice_layout,\n' +
    '           sp ∈ sp_elements где (sp.sp_type, sp.sp_element_id) =\n' +
    '                                  (layout.sp_type, layout.sp_element_id)}\n' +
    '       ∪ {копия строк осн. решётки с lattice_id ← sub.lattice_id и\n' +
    '          instance_id ≤ sub.max_gid\n' +
    '         | sub ∈ sub_lattices}'),

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
  body('Огибающая:'),
  code('A_i(t_k) = distanceToPolyline( x_i(t_k), design_poly )'),
  body('Линейная подгонка методом наименьших квадратов:'),
  code('A(t) ≈ m · t + A₀'),
  body('Экспоненциальная подгонка (по логарифму положительной части):'),
  code(
    'log A(t) ≈ γ · t + log A₀   ⇒   A(t) ≈ A₀ · exp(γ · t)\n' +
    'R²_exp = 1 − SS_res / SS_tot'),
  body('FFT при шаге выборки dt_sample = t_k − t_{k−1}:'),
  code(
    'Â(f) = FFT(A(t))\n' +
    'f_betatron = argmax_(f > f_min) |Â(f)|       (DC-бин маскируется)'),

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
