// Where the site and its repository live. GitHub Actions sets GITHUB_REPOSITORY, so a
// build in a fork links to the fork. A local build uses the default below.
const slug = process.env.GITHUB_REPOSITORY || 'dennisklein/FairRoot';
const owner = slug.split('/')[0];

/** The repository on GitHub, for example "https://github.com/dennisklein/FairRoot". */
export const repoUrl = `https://github.com/${slug}`;
/** The GitHub Pages host of the repository owner. */
export const site = `https://${owner.toLowerCase()}.github.io`;
/** The path of the site on the host. */
export const base = '/FairRoot';
/** The branch that the site documents and deploys from. */
export const branch = 'dev';
