// Одноразовый скрипт: пробегает docs/algorithms.md, добавляет \tag{X.Y}
// перед каждым закрывающим $$ (X — номер главы, Y — счётчик в главе).
// Запуск:  node docs/_tag_formulas.js
const fs = require('fs');
const path = require('path');
const file = path.join(__dirname, 'algorithms.md');
const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/);
let section = null, counter = 0, inFormula = false, tagged = 0;
const out = [];
for (const line of lines) {
  const m = line.match(/^## (\d+)\. /);
  if (m) { section = m[1]; counter = 0; out.push(line); continue; }
  if (line === '$$') {
    if (!inFormula) { inFormula = true; out.push(line); }
    else {
      if (section) {
        counter++; tagged++;
        out.push(String.raw`\tag{${section}.${counter}}`);
      }
      inFormula = false;
      out.push(line);
    }
    continue;
  }
  out.push(line);
}
fs.writeFileSync(file, out.join('\n'));
console.log('Inserted \\tag{} lines:', tagged);
