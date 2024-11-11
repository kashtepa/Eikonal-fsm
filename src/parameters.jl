# Общие для всех параметры
module Parameters

include("lense.jl")

using StaticArrays: SVector
using .Lenses

export lense_type
export R, n₀, center
export x₀
export T, t_start, t_stop
export h, A_tol, R_tol
export α_lim, β_lim
export x_lim, y_lim, z_lim
export φ_min, φ_max, ρ_min, ρ_max

"Радиус линзы"
const R = 1.0

"Коэффициент преломления среды"
const n₀ = 1.0

"Координаты центра линзы"
const center = SVector(0.0, 0.0, 0.0)

# Границы области для декартовых координат
const x_lim = (0.0, 5.0)
const y_lim = (-0.95, 0.95)
const z_lim = (-0.95, 0.95)

# Границы области для полярных координат
const φ_min =  0.0
const φ_max =  2π
const ρ_min = 0.0
const ρ_max =  2.0

# Диапазон направляющих углов
# в двумерном случае угол α фиксирован
const α_lim = Dict(
  "3d" => (0.1, π-0.1),
  "2d" => π/2,
)

const β_lim = Dict(
  "3d" => (0.1, pi-0.1),
  "2d" => (-π/2.01, π/2.01),
)


"Шаг для обычного метода Рунге-Кутты"
const h = 1.0e-2
"Точность для вложенного метода Рунге-Кутты"
const A_tol = 1.0e-7
const R_tol = 1.0e-7

"Временной отрезок"
const T = (start=0.0, stop=250h)

const t_start = T.start
const t_stop  = T.stop


"Начальные значения <=> позиция точечного источника"
const x₀ = SVector(-1.0, 0.0, 0)

# Начальные значения для полярных координат
const φ₀ = π
const ρ₀ = 1.0

end