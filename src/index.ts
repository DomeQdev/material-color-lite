const delinearized = (rgbComponent: number) => {
    const normalized = rgbComponent / 100.0;
    let delinearized = 0.0;
    if (normalized <= 0.0031308) {
        delinearized = normalized * 12.92;
    } else {
        delinearized = 1.055 * Math.pow(normalized, 1.0 / 2.4) - 0.055;
    }

    return Math.max(0, Math.min(255, Math.round(delinearized * 255.0)));
};

const linearized = (rgbComponent: number) => {
    const normalized = rgbComponent / 255.0;
    if (normalized <= 0.040449936) {
        return (normalized / 12.92) * 100.0;
    } else {
        return Math.pow((normalized + 0.055) / 1.055, 2.4) * 100.0;
    }
};

const argbFromRgb = (red: number, green: number, blue: number) => {
    return (0xff000000 | ((red & 255) << 16) | ((green & 255) << 8) | (blue & 255)) >>> 0;
};

const argbFromXyz = (x: number, y: number, z: number): number => {
    const r = delinearized(3.2413774792388685 * x + -1.5376652402851851 * y + -0.49885366846268053 * z);
    const g = delinearized(-0.9691452513005321 * x + 1.8758853451067872 * y + 0.04156585616912061 * z);
    const b = delinearized(0.05562093689691305 * x + -0.20395524564742123 * y + 1.0571799111220335 * z);

    return argbFromRgb(r, g, b);
};

const xyzFromArgb = (argb: number): number[] => {
    const r = linearized(redFromArgb(argb));
    const g = linearized(greenFromArgb(argb));
    const b = linearized(blueFromArgb(argb));

    return [
        0.41233895 * r + 0.35762064 * g + 0.18051042 * b,
        0.2126 * r + 0.7152 * g + 0.0722 * b,
        0.01932141 * r + 0.11916382 * g + 0.95034478 * b,
    ];
};

const e = 216.0 / 24389.0;
const kappa = 24389.0 / 27.0;

const labInvf = (ft: number) => {
    const ft3 = ft * ft * ft;

    return ft3 > e ? ft3 : (116 * ft - 16) / kappa;
};

const argbFromLstar = (lstar: number) => {
    const component = delinearized(yFromLstar(lstar));

    return argbFromRgb(component, component, component);
};

const labF = (t: number) => (t > e ? Math.cbrt(t) : (kappa * t + 16) / 116);
const lstarFromArgb = (argb: number) => 116.0 * labF(xyzFromArgb(argb)[1] / 100.0) - 16.0;
const yFromLstar = (lstar: number) => 100.0 * labInvf((lstar + 16.0) / 116.0);
const redFromArgb = (argb: number) => (argb >> 16) & 255;
const greenFromArgb = (argb: number) => (argb >> 8) & 255;
const blueFromArgb = (argb: number) => argb & 255;
const lerp = (start: number, stop: number, amount: number) => (1.0 - amount) * start + amount * stop;
const signum = (num: number): number => (num < 0 ? -1 : num > 0 ? 1 : 0);
const sanitizeDegreesDouble = (degrees: number) => (degrees + 360.0) % 360.0;

class ViewingConditions {
    static DEFAULT = ViewingConditions.make();

    static make(
        whitePoint = [95.047, 100.0, 108.883],
        adaptingLuminance = ((200.0 / Math.PI) * yFromLstar(50.0)) / 100.0,
        backgroundLstar = 50.0,
        surround = 2.0,
        discountingIlluminant = false
    ) {
        const rW = whitePoint[0] * 0.401288 + whitePoint[1] * 0.650173 + whitePoint[2] * -0.051461;
        const gW = whitePoint[0] * -0.250268 + whitePoint[1] * 1.204414 + whitePoint[2] * 0.045854;
        const bW = whitePoint[0] * -0.002079 + whitePoint[1] * 0.048952 + whitePoint[2] * 0.953127;
        const f = 0.8 + surround / 10.0;
        const c = f >= 0.9 ? lerp(0.59, 0.69, (f - 0.9) * 10.0) : lerp(0.525, 0.59, (f - 0.8) * 10.0);
        let d = discountingIlluminant
            ? 1.0
            : f * (1.0 - (1.0 / 3.6) * Math.exp((-adaptingLuminance - 42.0) / 92.0));
        d = d > 1.0 ? 1.0 : d < 0.0 ? 0.0 : d;
        const rgbD = [d * (100.0 / rW) + 1.0 - d, d * (100.0 / gW) + 1.0 - d, d * (100.0 / bW) + 1.0 - d];
        const k = 1.0 / (5.0 * adaptingLuminance + 1.0);
        const k4 = k * k * k * k;
        const k4F = 1.0 - k4;
        const fl = k4 * adaptingLuminance + 0.1 * k4F * k4F * Math.cbrt(5.0 * adaptingLuminance);
        const n = yFromLstar(backgroundLstar) / whitePoint[1];
        const z = 1.48 + Math.sqrt(n);
        const nbb = 0.725 / Math.pow(n, 0.2);
        const rgbAFactors = [
            Math.pow((fl * rgbD[0] * rW) / 100.0, 0.42),
            Math.pow((fl * rgbD[1] * gW) / 100.0, 0.42),
            Math.pow((fl * rgbD[2] * bW) / 100.0, 0.42),
        ];
        const rgbA = [
            (400.0 * rgbAFactors[0]) / (rgbAFactors[0] + 27.13),
            (400.0 * rgbAFactors[1]) / (rgbAFactors[1] + 27.13),
            (400.0 * rgbAFactors[2]) / (rgbAFactors[2] + 27.13),
        ];
        const aw = (2.0 * rgbA[0] + rgbA[1] + 0.05 * rgbA[2]) * nbb;

        return new ViewingConditions(n, aw, nbb, nbb, c, f, rgbD, fl, Math.pow(fl, 0.25), z);
    }

    private constructor(
        public n: number,
        public aw: number,
        public nbb: number,
        public ncb: number,
        public c: number,
        public nc: number,
        public rgbD: number[],
        public fl: number,
        public fLRoot: number,
        public z: number
    ) {}
}

class Cam16 {
    constructor(
        readonly hue: number,
        readonly chroma: number,
        readonly j: number,
        readonly q: number,
        readonly m: number,
        readonly s: number,
        readonly jstar: number,
        readonly astar: number,
        readonly bstar: number
    ) {}

    distance(other: Cam16): number {
        const dJ = this.jstar - other.jstar;
        const dA = this.astar - other.astar;
        const dB = this.bstar - other.bstar;

        return 1.41 * Math.pow(Math.sqrt(dJ * dJ + dA * dA + dB * dB), 0.63);
    }

    static fromInt(argb: number) {
        const viewingConditions = ViewingConditions.DEFAULT;
        const [x, y, z] = xyzFromArgb(argb);
        const rC = 0.401288 * x + 0.650173 * y - 0.051461 * z;
        const gC = -0.250268 * x + 1.204414 * y + 0.045854 * z;
        const bC = -0.002079 * x + 0.048952 * y + 0.953127 * z;
        const rD = viewingConditions.rgbD[0] * rC;
        const gD = viewingConditions.rgbD[1] * gC;
        const bD = viewingConditions.rgbD[2] * bC;
        const rAF = Math.pow((viewingConditions.fl * Math.abs(rD)) / 100.0, 0.42);
        const gAF = Math.pow((viewingConditions.fl * Math.abs(gD)) / 100.0, 0.42);
        const bAF = Math.pow((viewingConditions.fl * Math.abs(bD)) / 100.0, 0.42);
        const rA = (signum(rD) * 400.0 * rAF) / (rAF + 27.13);
        const gA = (signum(gD) * 400.0 * gAF) / (gAF + 27.13);
        const bA = (signum(bD) * 400.0 * bAF) / (bAF + 27.13);
        const a = (11.0 * rA - 12.0 * gA + bA) / 11.0;
        const bVal = (rA + gA - 2.0 * bA) / 9.0;
        const u = (20.0 * rA + 20.0 * gA + 21.0 * bA) / 20.0;
        const p2 = (40.0 * rA + 20.0 * gA + bA) / 20.0;
        const atanDegrees = (Math.atan2(bVal, a) * 180.0) / Math.PI;
        const hue = sanitizeDegreesDouble(atanDegrees);
        const hueRadians = (hue * Math.PI) / 180.0;
        const ac = p2 * viewingConditions.nbb;
        const j = 100.0 * Math.pow(ac / viewingConditions.aw, viewingConditions.c * viewingConditions.z);
        const q =
            (4.0 / viewingConditions.c) *
            Math.sqrt(j / 100.0) *
            (viewingConditions.aw + 4.0) *
            viewingConditions.fLRoot;
        const huePrime = hue < 20.14 ? hue + 360 : hue;
        const eHue = 0.25 * (Math.cos((huePrime * Math.PI) / 180.0 + 2.0) + 3.8);
        const p1 = (50000.0 / 13.0) * eHue * viewingConditions.nc * viewingConditions.ncb;
        const t = (p1 * Math.sqrt(a * a + bVal * bVal)) / (u + 0.305);
        const alpha = Math.pow(t, 0.9) * Math.pow(1.64 - Math.pow(0.29, viewingConditions.n), 0.73);
        const c = alpha * Math.sqrt(j / 100.0);
        const m = c * viewingConditions.fLRoot;
        const s = 50.0 * Math.sqrt((alpha * viewingConditions.c) / (viewingConditions.aw + 4.0));
        const jstar = ((1.0 + 100.0 * 0.007) * j) / (1.0 + 0.007 * j);
        const mstar = (1.0 / 0.0228) * Math.log(1.0 + 0.0228 * m);
        const astar = mstar * Math.cos(hueRadians);
        const bstar = mstar * Math.sin(hueRadians);

        return new Cam16(hue, c, j, q, m, s, jstar, astar, bstar);
    }

    static fromJch(j: number, c: number, h: number) {
        const viewingConditions = ViewingConditions.DEFAULT;
        const q =
            (4.0 / viewingConditions.c) *
            Math.sqrt(j / 100.0) *
            (viewingConditions.aw + 4.0) *
            viewingConditions.fLRoot;
        const m = c * viewingConditions.fLRoot;
        const alpha = c / Math.sqrt(j / 100.0);
        const s = 50.0 * Math.sqrt((alpha * viewingConditions.c) / (viewingConditions.aw + 4.0));
        const hueRadians = (h * Math.PI) / 180.0;
        const jstar = ((1.0 + 100.0 * 0.007) * j) / (1.0 + 0.007 * j);
        const mstar = (1.0 / 0.0228) * Math.log(1.0 + 0.0228 * m);
        const astar = mstar * Math.cos(hueRadians);
        const bstar = mstar * Math.sin(hueRadians);

        return new Cam16(h, c, j, q, m, s, jstar, astar, bstar);
    }

    toInt() {
        const viewingConditions = ViewingConditions.DEFAULT;
        const alpha = this.chroma === 0.0 || this.j === 0.0 ? 0.0 : this.chroma / Math.sqrt(this.j / 100.0);
        const t = Math.pow(alpha / Math.pow(1.64 - Math.pow(0.29, viewingConditions.n), 0.73), 1.0 / 0.9);
        const hRad = (this.hue * Math.PI) / 180.0;
        const eHue = 0.25 * (Math.cos(hRad + 2.0) + 3.8);
        const ac =
            viewingConditions.aw * Math.pow(this.j / 100.0, 1.0 / viewingConditions.c / viewingConditions.z);
        const p1 = eHue * (50000.0 / 13.0) * viewingConditions.nc * viewingConditions.ncb;
        const p2 = ac / viewingConditions.nbb;
        const hSin = Math.sin(hRad);
        const hCos = Math.cos(hRad);
        const gamma = (23.0 * (p2 + 0.305) * t) / (23.0 * p1 + 11.0 * t * hCos + 108.0 * t * hSin);
        const a = gamma * hCos;
        const b = gamma * hSin;
        const rA = (460.0 * p2 + 451.0 * a + 288.0 * b) / 1403.0;
        const gA = (460.0 * p2 - 891.0 * a - 261.0 * b) / 1403.0;
        const bA = (460.0 * p2 - 220.0 * a - 6300.0 * b) / 1403.0;
        const rCBase = Math.max(0, (27.13 * Math.abs(rA)) / (400.0 - Math.abs(rA)));
        const rC = signum(rA) * (100.0 / viewingConditions.fl) * Math.pow(rCBase, 1.0 / 0.42);
        const gCBase = Math.max(0, (27.13 * Math.abs(gA)) / (400.0 - Math.abs(gA)));
        const gC = signum(gA) * (100.0 / viewingConditions.fl) * Math.pow(gCBase, 1.0 / 0.42);
        const bCBase = Math.max(0, (27.13 * Math.abs(bA)) / (400.0 - Math.abs(bA)));
        const bC = signum(bA) * (100.0 / viewingConditions.fl) * Math.pow(bCBase, 1.0 / 0.42);
        const rF = rC / viewingConditions.rgbD[0];
        const gF = gC / viewingConditions.rgbD[1];
        const bF = bC / viewingConditions.rgbD[2];
        const x = 1.86206786 * rF - 1.01125463 * gF + 0.14918677 * bF;
        const y = 0.38752654 * rF + 0.62144744 * gF - 0.00897398 * bF;
        const z = -0.0158415 * rF - 0.03412294 * gF + 1.04996444 * bF;

        return argbFromXyz(x, y, z);
    }
}

const solveHctToInt = (hue: number, chroma: number, tone: number) => {
    if (chroma < 1.0 || Math.round(tone) <= 0.0 || Math.round(tone) >= 100.0) {
        return argbFromLstar(tone);
    }

    hue = sanitizeDegreesDouble(hue);

    let high = chroma;
    let mid = chroma;
    let low = 0.0;
    let isFirstLoop = true;

    let answer: Cam16 | null = null;

    while (Math.abs(low - high) >= 0.4) {
        const possibleAnswer = findCamByJ(hue, mid, tone);

        if (isFirstLoop) {
            if (possibleAnswer != null) {
                return possibleAnswer.toInt();
            } else {
                isFirstLoop = false;
                mid = low + (high - low) / 2.0;
                continue;
            }
        }

        if (possibleAnswer === null) {
            high = mid;
        } else {
            answer = possibleAnswer;
            low = mid;
        }

        mid = low + (high - low) / 2.0;
    }

    if (answer === null) {
        return argbFromLstar(tone);
    }
    return answer.toInt();
};

const findCamByJ = (hue: number, chroma: number, tone: number) => {
    let low = 0.0;
    let high = 100.0;
    let mid = 0.0;
    let bestdL = 1000.0;
    let bestdE = 1000.0;
    let bestCam: Cam16 | null = null;

    while (Math.abs(low - high) > 0.01) {
        mid = low + (high - low) / 2.0;
        const camBeforeGamut = Cam16.fromJch(mid, chroma, hue);
        const argb = camBeforeGamut.toInt();
        const camAfterGamut = Cam16.fromInt(argb);
        const Lstar = lstarFromArgb(argb);
        const dL = Math.abs(tone - Lstar);

        if (dL < 0.2) {
            const dE = camBeforeGamut.distance(camAfterGamut);
            if (dE <= 1.0) {
                bestdL = dL;
                bestdE = dE;
                bestCam = camAfterGamut;
            }
        }
        if (bestdL === 0.0 && bestdE === 0.0) {
            return bestCam;
        }
        if (Lstar < tone) {
            low = mid;
        } else {
            high = mid;
        }
    }

    return bestCam;
};

class Hct {
    private internalHue: number;
    private internalChroma: number;
    private argb: number;

    static from(hue: number, chroma: number, tone: number) {
        return new Hct(solveHctToInt(hue, chroma, tone));
    }

    static fromInt(argb: number) {
        return new Hct(argb);
    }

    toInt() {
        return this.argb;
    }

    get hue() {
        return this.internalHue;
    }

    get chroma() {
        return this.internalChroma;
    }

    private constructor(argb: number) {
        const { hue, chroma } = Cam16.fromInt(argb);
        this.internalHue = hue;
        this.internalChroma = chroma;
        this.argb = argb;
    }
}

class TonalPalette {
    private readonly cache = new Map<number, number>();

    static fromHueAndChroma(hue: number, chroma: number) {
        return new TonalPalette(hue, chroma);
    }

    private constructor(readonly hue: number, readonly chroma: number) {}

    tone(tone: number) {
        let argb = this.cache.get(tone);
        if (argb === undefined) {
            argb = Hct.from(this.hue, this.chroma, tone).toInt();
            this.cache.set(tone, argb);
        }
        return argb;
    }
}

class CorePalette {
    primary: TonalPalette;
    secondary: TonalPalette;
    tertiary: TonalPalette;
    neutral: TonalPalette;
    neutralVariant: TonalPalette;
    error: TonalPalette;

    static of(argb: number) {
        return new CorePalette(argb);
    }

    private constructor(argb: number) {
        const { hue, chroma } = Hct.fromInt(argb);

        this.primary = TonalPalette.fromHueAndChroma(hue, Math.max(48, chroma));
        this.secondary = TonalPalette.fromHueAndChroma(hue, 16);
        this.tertiary = TonalPalette.fromHueAndChroma(hue + 60, 24);
        this.neutral = TonalPalette.fromHueAndChroma(hue, 4);
        this.neutralVariant = TonalPalette.fromHueAndChroma(hue, 8);
        this.error = TonalPalette.fromHueAndChroma(25, 84);
    }
}

const argbToHex = (argb: number) => {
    return `#${(argb & 0x00ffffff).toString(16).padStart(6, "0")}`;
};

const argbToRgb = (argb: number) => {
    return `rgb(${redFromArgb(argb)}, ${greenFromArgb(argb)}, ${blueFromArgb(argb)})`;
};

const hexToArgb = (hex: string) => {
    const hexValue = hex.startsWith("#") ? hex.slice(1) : hex;

    const r = parseInt(hexValue.substring(0, 2), 16);
    const g = parseInt(hexValue.substring(2, 4), 16);
    const b = parseInt(hexValue.substring(4, 6), 16);

    return argbFromRgb(r, g, b);
};

export enum ColorRole {
    Primary = "primary",
    OnPrimary = "onPrimary",
    PrimaryContainer = "primaryContainer",
    OnPrimaryContainer = "onPrimaryContainer",
    Secondary = "secondary",
    OnSecondary = "onSecondary",
    SecondaryContainer = "secondaryContainer",
    OnSecondaryContainer = "onSecondaryContainer",
    Tertiary = "tertiary",
    OnTertiary = "onTertiary",
    TertiaryContainer = "tertiaryContainer",
    OnTertiaryContainer = "onTertiaryContainer",
    Error = "error",
    OnError = "onError",
    ErrorContainer = "errorContainer",
    OnErrorContainer = "onErrorContainer",
    Background = "background",
    OnBackground = "onBackground",
    Surface = "surface",
    OnSurface = "onSurface",
    SurfaceVariant = "surfaceVariant",
    OnSurfaceVariant = "onSurfaceVariant",
    Outline = "outline",
    OutlineVariant = "outlineVariant",
    Shadow = "shadow",
    Scrim = "scrim",
    InverseSurface = "inverseSurface",
    InverseOnSurface = "inverseOnSurface",
    InversePrimary = "inversePrimary",
}

type SchemeMap = Record<ColorRole, (core: CorePalette) => number>;

const lightSchemeMap: SchemeMap = {
    [ColorRole.Primary]: (core) => core.primary.tone(40),
    [ColorRole.OnPrimary]: (core) => core.primary.tone(100),
    [ColorRole.PrimaryContainer]: (core) => core.primary.tone(90),
    [ColorRole.OnPrimaryContainer]: (core) => core.primary.tone(10),
    [ColorRole.Secondary]: (core) => core.secondary.tone(40),
    [ColorRole.OnSecondary]: (core) => core.secondary.tone(100),
    [ColorRole.SecondaryContainer]: (core) => core.secondary.tone(90),
    [ColorRole.OnSecondaryContainer]: (core) => core.secondary.tone(10),
    [ColorRole.Tertiary]: (core) => core.tertiary.tone(40),
    [ColorRole.OnTertiary]: (core) => core.tertiary.tone(100),
    [ColorRole.TertiaryContainer]: (core) => core.tertiary.tone(90),
    [ColorRole.OnTertiaryContainer]: (core) => core.tertiary.tone(10),
    [ColorRole.Error]: (core) => core.error.tone(40),
    [ColorRole.OnError]: (core) => core.error.tone(100),
    [ColorRole.ErrorContainer]: (core) => core.error.tone(90),
    [ColorRole.OnErrorContainer]: (core) => core.error.tone(10),
    [ColorRole.Background]: (core) => core.neutral.tone(99),
    [ColorRole.OnBackground]: (core) => core.neutral.tone(10),
    [ColorRole.Surface]: (core) => core.neutral.tone(99),
    [ColorRole.OnSurface]: (core) => core.neutral.tone(10),
    [ColorRole.SurfaceVariant]: (core) => core.neutralVariant.tone(90),
    [ColorRole.OnSurfaceVariant]: (core) => core.neutralVariant.tone(30),
    [ColorRole.Outline]: (core) => core.neutralVariant.tone(50),
    [ColorRole.OutlineVariant]: (core) => core.neutralVariant.tone(80),
    [ColorRole.Shadow]: (core) => core.neutral.tone(0),
    [ColorRole.Scrim]: (core) => core.neutral.tone(0),
    [ColorRole.InverseSurface]: (core) => core.neutral.tone(20),
    [ColorRole.InverseOnSurface]: (core) => core.neutral.tone(95),
    [ColorRole.InversePrimary]: (core) => core.primary.tone(80),
};

const darkSchemeMap: SchemeMap = {
    [ColorRole.Primary]: (core) => core.primary.tone(80),
    [ColorRole.OnPrimary]: (core) => core.primary.tone(20),
    [ColorRole.PrimaryContainer]: (core) => core.primary.tone(30),
    [ColorRole.OnPrimaryContainer]: (core) => core.primary.tone(90),
    [ColorRole.Secondary]: (core) => core.secondary.tone(80),
    [ColorRole.OnSecondary]: (core) => core.secondary.tone(20),
    [ColorRole.SecondaryContainer]: (core) => core.secondary.tone(30),
    [ColorRole.OnSecondaryContainer]: (core) => core.secondary.tone(90),
    [ColorRole.Tertiary]: (core) => core.tertiary.tone(80),
    [ColorRole.OnTertiary]: (core) => core.tertiary.tone(20),
    [ColorRole.TertiaryContainer]: (core) => core.tertiary.tone(30),
    [ColorRole.OnTertiaryContainer]: (core) => core.tertiary.tone(90),
    [ColorRole.Error]: (core) => core.error.tone(80),
    [ColorRole.OnError]: (core) => core.error.tone(20),
    [ColorRole.ErrorContainer]: (core) => core.error.tone(30),
    [ColorRole.OnErrorContainer]: (core) => core.error.tone(80),
    [ColorRole.Background]: (core) => core.neutral.tone(10),
    [ColorRole.OnBackground]: (core) => core.neutral.tone(90),
    [ColorRole.Surface]: (core) => core.neutral.tone(10),
    [ColorRole.OnSurface]: (core) => core.neutral.tone(90),
    [ColorRole.SurfaceVariant]: (core) => core.neutralVariant.tone(30),
    [ColorRole.OnSurfaceVariant]: (core) => core.neutralVariant.tone(80),
    [ColorRole.Outline]: (core) => core.neutralVariant.tone(60),
    [ColorRole.OutlineVariant]: (core) => core.neutralVariant.tone(30),
    [ColorRole.Shadow]: (core) => core.neutral.tone(0),
    [ColorRole.Scrim]: (core) => core.neutral.tone(0),
    [ColorRole.InverseSurface]: (core) => core.neutral.tone(90),
    [ColorRole.InverseOnSurface]: (core) => core.neutral.tone(20),
    [ColorRole.InversePrimary]: (core) => core.primary.tone(40),
};
abstract class SchemeBase {
    private readonly corePalette: CorePalette;
    private readonly returnRGB: boolean;
    protected abstract readonly map: SchemeMap;

    constructor(sourceHexColor: string, returnRGB = false) {
        this.corePalette = CorePalette.of(hexToArgb(sourceHexColor));
        this.returnRGB = returnRGB;
    }

    protected _get(role: ColorRole): string {
        const color = this.map[role](this.corePalette);

        return this.returnRGB ? argbToRgb(color) : argbToHex(color);
    }

    getTone(tonalPalette: keyof CorePalette, tone: number): string {
        const color = this.corePalette[tonalPalette].tone(tone);

        return this.returnRGB ? argbToRgb(color) : argbToHex(color);
    }

    get primary() {
        return this._get(ColorRole.Primary);
    }
    get onPrimary() {
        return this._get(ColorRole.OnPrimary);
    }
    get primaryContainer() {
        return this._get(ColorRole.PrimaryContainer);
    }
    get onPrimaryContainer() {
        return this._get(ColorRole.OnPrimaryContainer);
    }
    get secondary() {
        return this._get(ColorRole.Secondary);
    }
    get onSecondary() {
        return this._get(ColorRole.OnSecondary);
    }
    get secondaryContainer() {
        return this._get(ColorRole.SecondaryContainer);
    }
    get onSecondaryContainer() {
        return this._get(ColorRole.OnSecondaryContainer);
    }
    get tertiary() {
        return this._get(ColorRole.Tertiary);
    }
    get onTertiary() {
        return this._get(ColorRole.OnTertiary);
    }
    get tertiaryContainer() {
        return this._get(ColorRole.TertiaryContainer);
    }
    get onTertiaryContainer() {
        return this._get(ColorRole.OnTertiaryContainer);
    }
    get error() {
        return this._get(ColorRole.Error);
    }
    get onError() {
        return this._get(ColorRole.OnError);
    }
    get errorContainer() {
        return this._get(ColorRole.ErrorContainer);
    }
    get onErrorContainer() {
        return this._get(ColorRole.OnErrorContainer);
    }
    get background() {
        return this._get(ColorRole.Background);
    }
    get onBackground() {
        return this._get(ColorRole.OnBackground);
    }
    get surface() {
        return this._get(ColorRole.Surface);
    }
    get onSurface() {
        return this._get(ColorRole.OnSurface);
    }
    get surfaceVariant() {
        return this._get(ColorRole.SurfaceVariant);
    }
    get onSurfaceVariant() {
        return this._get(ColorRole.OnSurfaceVariant);
    }
    get outline() {
        return this._get(ColorRole.Outline);
    }
    get outlineVariant() {
        return this._get(ColorRole.OutlineVariant);
    }
    get shadow() {
        return this._get(ColorRole.Shadow);
    }
    get scrim() {
        return this._get(ColorRole.Scrim);
    }
    get inverseSurface() {
        return this._get(ColorRole.InverseSurface);
    }
    get inverseOnSurface() {
        return this._get(ColorRole.InverseOnSurface);
    }
    get inversePrimary() {
        return this._get(ColorRole.InversePrimary);
    }
}

export class LightScheme extends SchemeBase {
    protected readonly map = lightSchemeMap;
}

export class DarkScheme extends SchemeBase {
    protected readonly map = darkSchemeMap;
}
