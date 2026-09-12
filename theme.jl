# Shared figure style. Every script in this repository includes this file, so
# that figures across the whole archive share fonts, sizing and axis treatment.

using CairoMakie, MathTeXEngine

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
