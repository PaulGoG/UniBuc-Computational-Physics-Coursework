# Shared figure style. Every script in this repository includes this file, so
# that figures across the whole archive share fonts, sizing, axis treatment and
# the small set of formatting helpers below.

using CairoMakie, MathTeXEngine, LaTeXStrings

set_theme!(Theme(
    fonts = (;
        regular = texfont(:text),
        bold = texfont(:bold),
        italic = texfont(:italic),
    ),
    fontsize = 22,
    figure_padding = 16,
    Axis = (
        xgridstyle = :dash, ygridstyle = :dash,
        xgridcolor = (:grey, 0.12), ygridcolor = (:grey, 0.12),
        xminorticksvisible = false, yminorticksvisible = false,
        xtickalign = 1, ytickalign = 1,
    ),
))

"Okabe–Ito, colourblind-safe. One consistent colour per quantity across the repository."
const PALETTE = (
    blue   = "#0072B2",
    orange = "#E69F00",
    green  = "#009E73",
    red    = "#D55E00",
    purple = "#CC79A7",
    sky    = "#56B4E9",
    yellow = "#F0E442",
    black  = "#000000",
)

"""
Marker sizes, one scale for the whole repository, chosen against the shared
`fontsize = 22` so that symbols stay legible when a figure is reduced to
journal column width.

  * `cloud` — scatter of thousands of points, where larger markers would only
    deepen the overplotting
  * `dense` — series with enough points that full-size markers would collide
  * `data` — the default for a measured or computed series
  * `emphasis` — single highlighted points: a fixed point, a fitted optimum
  * `key` — markers built by hand for a legend, which must stay readable even
    where the series itself is drawn small
"""
const MARKERSIZE = (cloud = 5, dense = 10, data = 14, emphasis = 19, key = 16)

_decimal_label(e::Integer) =
    e >= 0 ? latexstring(string(10^e)) : latexstring("0." * "0"^(-e - 1) * "1")

_power_label(e::Integer) = e == 0 ? L"1" : e == 1 ? L"10" : latexstring("10^{$e}")

"""
    logticks(e_min, e_max; step = 1, style = :auto)

Explicit decade ticks for a base-10 logarithmic axis running from `10^e_min` to
`10^e_max`, returned as the `(values, labels)` tuple an `Axis` takes.

Spans of four decades or fewer are labelled as plain decimals (`0.01, 0.1, 1,
10`) provided the decades themselves stay within `10^-4` to `10^4`, beyond which
a decimal is a run of zeros and the power of ten is the readable form; wider
spans use powers of ten, with `10^0` written `1` and `10^1` written `10`. Pass
`style = :decimal` or `style = :power` to override the choice.

Every logarithmic axis in this repository sets its ticks through this function.
Left to itself Makie labels the unit decade as `10^0` on some axes and as `1` on
others, and on ranges narrower than a decade it falls back to fractional
exponents such as `10^{-0.4}`, which is not a form any of these quantities are
read in.

# Example

```julia
Axis(fig[1, 1]; xscale = log10, xticks = logticks(-4, 0))
```
"""
function logticks(e_min::Integer, e_max::Integer; step::Integer = 1, style::Symbol = :auto)
    e_min <= e_max || throw(ArgumentError("e_min = $e_min exceeds e_max = $e_max"))
    step > 0 || throw(ArgumentError("step must be positive, got $step"))
    style in (:auto, :decimal, :power) ||
        throw(ArgumentError("style must be :auto, :decimal or :power, got :$style"))
    exponents = collect(e_min:step:e_max)
    readable = e_max - e_min <= 4 && e_min >= -4 && e_max <= 4
    decimal = style == :decimal || (style == :auto && readable)
    labels = [decimal ? _decimal_label(e) : _power_label(e) for e in exponents]
    return ([10.0^e for e in exponents], labels)
end

"""
    it(s)

`s` set in the italic of the shared Computer Modern family, for use inside a
`rich(...)` run: `rich("Lorenz: ", it("σ"), " = 10")`.

`rich` rather than a `LaTeXString` wherever a label carries ordinary prose,
because MathTeXEngine parses everything in an `L"…"` — including the text
outside `\$…\$` — as mathematics, so a colon picks up relation spacing and
renders as `Lorenz : σ = 10`. Reserve `L"…"` for labels that are entirely
mathematical.
"""
it(s) = rich(s, font = :italic)

_trim(v::Real) = (s = string(v); endswith(s, ".0") ? s[1:(end - 2)] : s)

"""
    sci(x; digits = 2)

`x` written as a LaTeX power of ten — `sci(1.8e-7)` gives `1.8 \\times 10^{-7}`.

Exponents of −1, 0 and 1 are folded into the mantissa (`0.18`, `1.8`, `18`) and a
mantissa of 1 is dropped, so the result never reads `1 \\times 10^{n}` or carries
a `\\times 10^{0}`. Figures and papers never show computer notation such as
`1.8e-07`, which is what `@sprintf("%e", …)` produces.
"""
function sci(x::Real; digits::Integer = 2)
    iszero(x) && return L"0"
    isfinite(x) || throw(ArgumentError("sci expects a finite number, got $x"))
    e = floor(Int, log10(abs(x)))
    m = round(x / 10.0^e; digits = digits)
    if abs(m) >= 10                      # rounding can carry the mantissa to 10
        m /= 10
        e += 1
    end
    -1 <= e <= 1 && return latexstring(_trim(round(x; sigdigits = digits + 1)))
    mantissa = _trim(m)
    mantissa == "1" && return latexstring("10^{$e}")
    mantissa == "-1" && return latexstring("-10^{$e}")
    return latexstring("$mantissa \\times 10^{$e}")
end

"""
    savefigure(fig, dir, name)

Write `fig` to `dir/name.png` well above 300 dpi, creating `dir` if needed.
Callers pass `joinpath(@__DIR__, "figures")` so that nothing depends on the
working directory. Returns the path written.
"""
function savefigure(fig, dir::AbstractString, name::AbstractString)
    mkpath(dir)
    path = joinpath(dir, name * ".png")
    save(path, fig; px_per_unit = 4)
    return path
end
