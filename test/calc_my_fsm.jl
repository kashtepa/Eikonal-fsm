include("../src/parameters.jl")
include("../src/lense.jl")
include("../fsm/fsm.jl")

using StaticArrays: SVector

using .FSM: fsm
using .Lenses

# Используем собственную реализацию FSM

const h = 0.01
const NSweeps = 10

const G = zeros(Bool, 1000, 1000)

# Центр линзы
const centre = (5.0, 5.0)
const R = 3.0
# Положение точечного источника
const source = (centre .- R / sqrt(2))
const LENSE = Maxwell(R, 1.0, SVector(centre..., 0.0))
# Точечный источник, на диагонали
# чтобы получить индексы мы делим на шаг h
G[round.(Int, (centre .- R / sqrt(2)) ./ h)...] = true

function nf(x, y)
  return Lenses.n(SVector(x, y, 0.0), LENSE)
end

X, Y, u = fsm(nf, G, h, NSweeps)


using CairoMakie

fig01 = Figure(size=(500, 500))
ax01 = Axis(fig01[1, 1])

# Рисуем контур линзы виде окружности
const lense_contour = Circle(Point2(centre...), R)
lines!(ax01, lense_contour, color=:blue)

contour!(ax01, X, Y, u, levels=100)
scatter!(ax01, source..., color=:red)

save("fsm.png", fig01)

#=
Получаемая картинка совпадает с библиотекой fsm
но происходит дублирование источников, с чем надо разобраться.
Также возникло сомнение в правильности рисования фронтов 
при решении методом характеристик
=#