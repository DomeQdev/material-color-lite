# 🚀 Material Color Lite ✨

A lightweight, high-performance, and tree-shakable alternative to `@material/material-color-utilities` for generating Material Design 3 color schemes.

### Why Material Color Lite?

*   **⚡ Blazing Fast:** Optimized for performance with zero dependencies.
*   **📦 Tiny & Tree-Shakable:** Your bundle includes only the code you use.
*   **🎨 Simple API:** Generate exactly the color roles you need with a clean, functional API.

---

## Quickstart

### 1. Installation

```bash
npm install material-color-lite
```

### 2. Usage

Generate light and dark theme colors from a single source color. The library returns only the roles you request, keeping your code clean and efficient.

```typescript
import { generateDarkScheme, generateLightScheme, ColorRole } from "material-color-lite";

const sourceColor = "#6750A4";

const requiredRoles = [ // Define which color roles you need for your theme
    ColorRole.Primary,
    ColorRole.OnPrimary,
    ColorRole.PrimaryContainer,
    ColorRole.Surface,
    ColorRole.OnSurface,
    ColorRole.Error,
];

const lightTheme = generateLightScheme(sourceColor, requiredRoles);
/*
{
  primary: '#7e67bd',
  onPrimary: '#ffffff',
  primaryContainer: '#f3d9ff',
  surface: '#fffbff',
  onSurface: '#212023',
  error: '#c04e44'
}
*/

const darkTheme = generateDarkScheme(sourceColor, requiredRoles);
/*
{
  primary: '#ddc4ff',
  onPrimary: '#462e81',
  primaryContainer: '#634ca0',
  surface: '#212023',
  onSurface: '#ede9ea',
  error: '#ff8b7a'
}
*/
```
