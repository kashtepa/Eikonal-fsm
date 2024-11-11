include("../src/lense.jl")

# Используем метод FSM, реализованный в библиотеке Eikonal

using Eikonal
using StaticArrays: SVector
using .Lenses

# Параметрическое уравнение окружности
circle_xy(center, R, φ) = round.(Int, center .+ (R*cos(φ), R*sin(φ)))

# Все координаты задаем в виде целых чисел, так как используется целочисленная сетка

# Центр линзы
const center = (500, 500)
# Радиус линзы
const R = 300
# Коэффициент преломления среды
const n₀ = 2.0
# Позиция источника относительно линзы
const source = circle_xy(center, R, pi)

const LENSE = select_lense(length(ARGS)>=1 ? ARGS[1] : "")(R, n₀, SVector(center..., 0.0))
const RAYS = false

println("Выбрана линза $(LENSE.name)")

# Размер сетки
const tsize = (1000, 1000)

fsm = FastSweeping(Float64, tsize)

# Вычисление коэффициента преломления в точках сетки
for x=1:tsize[1], y=1:tsize[2]
    fsm.v[x, y] = Lenses.n(SVector(x, y, 0.0), LENSE)
end

init!(fsm, source)

println("Подметание")
sweep!(fsm, verbose=false)

println("Рисование")
using CairoMakie

fig01 = Figure(fontsize=18, pt_per_unit=1)
ax01 = Axis(fig01[1, 1], xlabel = L"$x_n$", ylabel = L"$y_n$")

ax01.aspect = DataAspect()
ax01.autolimitaspect = 1
# colsize!(fig01.layout, 1, Aspect(1, 1))
rowsize!(fig01.layout, 1, Aspect(1, 1))

# Контур линзы виде окружности
const lense_contour = Circle(Point2(center .|> Float64), R)

# Рисуем контур линзы
lines!(ax01, lense_contour, color=:blue)

# Фронты в виде контуров функции от двух переменных u(x, y)
contour!(ax01, fsm.t, levels=100, colormap=:grays)

# Источник в виде точки
scatter!(ax01, source..., color=:red)

# Вычисление лучей
if RAYS
  #for θ = [-pi/6, -pi/24, -pi/48, -pi/960, 0, pi/960, pi/300, pi/200, pi/48, pi/24, pi/6]
  # pos = circle_xy(center, R, θ)
  for pos in ((900, 300), (900, 400), (900, 500),
    (900, 600), (900, 700), (900, 800))
    r1 = ray(fsm.t, pos)
    lines!(ax01, r1, color=:green)
  end
end


println("Сохранение")
if RAYS
  save("img/FSM/Линза_$(LENSE.name)_фронты_и_лучи.pdf", fig01)
else
  save("img/FSM/Линза_$(LENSE.name)_фронты.pdf", fig01)
end
