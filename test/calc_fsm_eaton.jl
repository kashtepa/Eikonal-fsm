include("../src/lense.jl")

# Используем метод FSM, реализованный в библиотеке Eikonal
# для расчета линзы Итона. Местоположение источников и сама
# конфигурация линзы для Итона настолько отличается, что 
# имеет смысл написать отдельную программу

using Eikonal
using StaticArrays: SVector
using .Lenses

# Параметрическое уравнение окружности
circle_xy(center, R, φ) = round.(Int, center .+ (R*cos(φ), R*sin(φ)))

# Все координаты задаем в виде целых чисел, так как используется целочисленная сетка

# Размер сетки
const tsize = (1000, 1000)
# Центр линзы
const center = (500, 500)
# Радиус линзы
const R = 250
# Коэффициент преломления среды
const n₀ = 1.0
# Позиция источника относительно центра линзы
const source = (center[1]+125, center[2])
# Местонахождение фокуса (симметрично относительно центра линзы)
const focus =  (center[1]-125, center[2])

const LENSE = Eaton(R, n₀, SVector(center..., 0.0))
const RAYS = true

println("Выбрана линза $(LENSE.name)")


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
const lense_contour01 = Circle(Point2(center .|> Float64), R)
const lense_contour02 = Circle(Point2(center .|> Float64), 2R)

# Рисуем контур линзы
lines!(ax01, lense_contour01, color=:blue)
lines!(ax01, lense_contour02, color=:blue)

# Фронты в виде контуров функции от двух переменных u(x, y)
contour!(ax01, fsm.t, levels=100, colormap=:grays)

# Для линзы Итона адекватно восстановить лучи не получается
if RAYS
  #for θ = [-pi/6, -pi/24, -pi/48, -pi/960, 0, pi/960, pi/300, pi/200, pi/48, pi/24, pi/6]
  # pos = circle_xy(center, R, θ)
  for pos in ((400, 300), (400, 400), (400, 500),
    (400, 600), (400, 700), (400, 800))
    r1 = ray(fsm.t, pos)
    lines!(ax01, r1, color=:green)
  end
end

# Источники в виде точек
scatter!(ax01, source..., color=:red)
scatter!(ax01, focus..., color=:blue)

println("Сохранение")
if RAYS
  save("img/FSM/Линза_$(LENSE.name)_фронты_и_лучи.pdf", fig01)
else
  save("img/FSM/Линза_$(LENSE.name)_фронты.pdf", fig01)
end
