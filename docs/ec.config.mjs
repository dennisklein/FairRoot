// Expressive Code options for the code blocks and the <Code> component.
// They are in their own file because the line-number plugin is not serializable to JSON.
import { defineEcConfig } from '@astrojs/starlight/expressive-code';
import { pluginLineNumbers } from '@expressive-code/plugin-line-numbers';
import codeTheme from './src/styles/code-theme.json' with { type: 'json' };

export default defineEcConfig({
  themes: [codeTheme],
  // Line numbers only where a block asks for them (the code issue pages).
  plugins: [pluginLineNumbers()],
  defaultProps: { showLineNumbers: false },
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
});
