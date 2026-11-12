// The register of known problems in the FairRoot code (src/data/code-issues.json).
import data from '../data/code-issues.json';
import { base } from './model.mjs';
import { SOURCE_DIRS, SOURCE_EXT } from '../../scripts/api.mjs';

export const SEVERITIES = [
  { id: 'high', label: 'High', text: 'A crash, memory corruption, data loss or a wrong physics result in a common setup.' },
  { id: 'medium', label: 'Medium', text: 'A wrong result, leak or trap in a less common setup, or a build problem.' },
  { id: 'low', label: 'Low', text: 'A wrong message, comment or test, dead code, or a small inconsistency.' },
];

export const CATEGORIES = {
  crash: 'Crash',
  'wrong-result': 'Wrong result',
  memory: 'Memory',
  'api-trap': 'API trap',
  build: 'Build',
  message: 'Message',
  'dead-code': 'Dead code',
  comment: 'Comment',
  test: 'Test',
  docs: 'Docs',
};

export const commit = data.commit;
export const issues = data.issues;
export const issueUrl = (id) => `${base}/issues/${id.toLowerCase()}/`;

/** The area of a file, for the filter: fairroot/<library>, examples/<group>, templates, cmake, … */
export function areaOf(file) {
  const p = file.split('/');
  if (p[0] === 'fairroot' && p.length > 2) return `fairroot/${p[1]}`;
  if (p[0] === 'examples' && p.length > 2) return `examples/${p[1]}`;
  if (p.length === 1) return '(top level)';
  return p[0];
}

const inBrowser = (file) => SOURCE_DIRS.includes(file.split('/')[0]) && SOURCE_EXT.test(file);
const range = (i) => (i.line ? `#L${i.line}${i.lineEnd ? `-L${i.lineEnd}` : ''}` : '');

/** Link to the file in the source browser (current dev), or to GitHub for other files. */
export function codeLink(i, repo) {
  return inBrowser(i.file) ? `${base}/source/${i.file}/${range(i)}` : `${repo}/blob/dev/${i.file}${range(i)}`;
}
/** Link to the lines on GitHub at the commit where the problem was found. */
export const foundLink = (i, repo) => `${repo}/blob/${commit}/${i.file}${range(i)}`;
export const where = (i) => `${i.file}${i.line ? `:${i.line}${i.lineEnd ? `–${i.lineEnd}` : ''}` : ''}`;
