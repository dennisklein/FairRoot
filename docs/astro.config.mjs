// @ts-check
import { readFileSync } from 'node:fs';
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import starlightLinksValidator from 'starlight-links-validator';
import { satteri } from '@astrojs/markdown-satteri';
import baseLinks from './scripts/base-links.mjs';
import nowrapCode from './scripts/nowrap-code.mjs';
import fairrootData from './scripts/integration.mjs';
import { site, base, repoUrl, branch } from './scripts/site.mjs';

const codeTheme = JSON.parse(readFileSync(new URL('./src/styles/code-theme.json', import.meta.url), 'utf8'));

export default defineConfig({
  site,
  base,
  markdown: {
    processor: satteri({ hastPlugins: [baseLinks({ base }), nowrapCode()] }),
    // Keep -- and ... as typed: options like -DBUILD_EXAMPLES and quoted log output must stay exact.
    smartypants: false,
  },
  integrations: [
    fairrootData(),
    starlight({
      title: 'FairRoot',
      description:
        'FairRoot is a C++ framework for simulation, reconstruction and analysis of particle physics experiments. It is based on ROOT and the Virtual Monte Carlo.',
      favicon: '/favicon.svg',
      social: [{ icon: 'github', label: 'GitHub', href: repoUrl }],
      editLink: { baseUrl: `${repoUrl}/edit/${branch}/docs/` },
      lastUpdated: true,
      disable404Route: true,
      customCss: [
        '@fontsource-variable/archivo/wdth.css',
        '@fontsource-variable/jetbrains-mono/wght.css',
        '@fontsource-variable/jetbrains-mono/wght-italic.css',
        './src/styles/theme.css',
      ],
      components: {
        Header: './src/components/overrides/Header.astro',
        PageTitle: './src/components/overrides/PageTitle.astro',
        Footer: './src/components/overrides/Footer.astro',
        PageSidebar: './src/components/overrides/PageSidebar.astro',
        Hero: './src/components/overrides/Hero.astro',
        MobileMenuFooter: './src/components/overrides/MobileMenuFooter.astro',
      },
      expressiveCode: {
        themes: [codeTheme],
        useStarlightDarkModeSwitch: false,
        useStarlightUiThemeColors: false,
        styleOverrides: {
          borderRadius: '6px',
          borderColor: 'var(--fr-code-rule)',
          codeBackground: 'var(--fr-code-bg)',
          codeFontFamily: 'var(--fr-font-mono)',
          codeFontSize: '0.84rem',
          codeLineHeight: '1.65',
          uiFontFamily: 'var(--fr-font-body)',
          frames: {
            editorBackground: 'var(--fr-code-bg)',
            editorTabBarBackground: 'var(--fr-code-bg)',
            editorActiveTabBackground: 'var(--fr-code-bg)',
            editorActiveTabIndicatorTopColor: 'transparent',
            editorActiveTabIndicatorBottomColor: 'var(--fr-saffron)',
            editorTabBarBorderBottomColor: 'var(--fr-code-rule)',
            terminalBackground: 'var(--fr-code-bg)',
            terminalTitlebarBackground: 'var(--fr-code-bg)',
            terminalTitlebarBorderBottomColor: 'var(--fr-code-rule)',
            terminalTitlebarDotsForeground: 'var(--fr-stage-line)',
            terminalTitlebarDotsOpacity: '1',
            frameBoxShadowCssValue: 'none',
            inlineButtonBorder: 'var(--fr-code-rule)',
          },
        },
      },
      sidebar: [
        { label: 'Start', items: [{ autogenerate: { directory: 'start' } }] },
        { label: 'Concepts', items: [{ autogenerate: { directory: 'concepts' } }] },
        { label: 'Simulation', items: [{ autogenerate: { directory: 'simulation' } }] },
        { label: 'Reconstruction', items: [{ autogenerate: { directory: 'reconstruction' } }] },
        { label: 'Online and FairMQ', collapsed: true, items: [{ autogenerate: { directory: 'online' } }] },
        { label: 'Tools', collapsed: true, items: [{ autogenerate: { directory: 'tools' } }] },
        { label: 'Reference', items: [{ autogenerate: { directory: 'reference' } }] },
        { label: 'Examples', collapsed: true, items: [{ autogenerate: { directory: 'examples' } }] },
        { label: 'Project', collapsed: true, items: [{ autogenerate: { directory: 'project' } }] },
        {
          label: 'Code',
          items: [
            { label: 'API reference', link: '/api/' },
            { label: 'Source browser', link: '/source/' },
            { label: 'Code issues', link: '/issues/' },
          ],
        },
      ],
      plugins: [
        starlightLinksValidator({
          errorOnRelativeLinks: false,
          // Generated pages. scripts/check-links.mjs checks these links in dist/.
          // Links to /repo/<path> go to GitHub (scripts/base-links.mjs).
          exclude: [`${base}/api/**`, `${base}/source/**`, `${base}/issues/**`, '/repo/**', `${base}/repo/**`],
        }),
      ],
    }),
  ],
});
