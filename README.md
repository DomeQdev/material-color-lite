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

Generate light and dark theme colors from a single source color. The library uses a lazy-loading approach, so colors are only generated when you access them.

```typescript
import { DarkScheme, LightScheme } from "material-color-lite";

const sourceColor = "#6750A4";

const lightTheme = new LightScheme(sourceColor);
const darkTheme = new DarkScheme(sourceColor);

console.log(lightTheme.primary); // '#7e67bd'
console.log(darkTheme.primary); // '#ddc4ff'

console.log(lightTheme.surface); // '#fffbff'
console.log(darkTheme.surface); // '#212023'
```
