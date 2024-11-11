module Lenses
  # using Base: explicit_manifest_deps_get
  using StaticArrays: SVector
  using LinearAlgebra: norm

  export Maxwell, Luneburg, Eaton
  export select_lense
  export eikonal_equation, geodesic_equation
  export refractive_index_and_its_derivative

  "Просто какая-то линза"
  abstract type Lense end

  "Линза Максвелла"
  struct Maxwell{T<:Real, V<:Real} <: Lense
    "Название линзы"
    name::String
    "Радиус линзы"
    R::T
    "Коэффициент преломления среды"
    n₀::T
    "Центр линзы"
    center::SVector{3, V}
  end

  "Линза Люниберга"
  struct Luneburg{T<:Real, V<:Real} <: Lense
    "Название линзы"
    name::String
    "Радиус линзы"
    R::T
    "Коэффициент преломления среды"
    n₀::T
    "Центр линзы"
    center::SVector{3, V}
  end

  "Линза Итона"
  struct Eaton{T<:Real, V<:Real} <: Lense
    "Название линзы"
    name::String
    "Радиус линзы"
    R::T
    "Коэффициент преломления среды"
    n₀::T
    "Центр линзы"
    center::SVector{3, V}
  end

  function Maxwell(R, n₀, center)
    Maxwell("Максвелла", promote(R, n₀)..., SVector(center))
  end

  function Luneburg(R, n₀, center)
    Luneburg("Люниберга", promote(R, n₀)..., SVector(center))
  end

  function Eaton(R, n₀, center)
    Eaton("Итона", promote(R, n₀)..., SVector(center))
  end


  "Выбрать поддерживаемую линзу по имени, по умолчанию Maxwell"
  function select_lense(name::String)
    if name == "maxwell"
      return Maxwell
    elseif name == "luneburg"
      return Luneburg
    elseif name == "eaton"
      return Eaton
    else
      return Maxwell
    end
  end

  "Расстояние от центра линзы до точки в пространстве"
  r(x::SVector{3, T}, L::Lense) where {T<:Real} = norm(x - L.center)

  """Коэффициент преломления линзы Люниберга"""
  function n(x::SVector{3, T}, L::Luneburg)::Float64 where {T<:Real}
    if r(x, L) <= L.R
      return L.n₀ * sqrt(2 - (r(x, L)/L.R)^2)
    else
      return L.n₀
    end
  end


  """Частные производные коэффициента преломления линзы Люниберга"""
  function ∂n(x::SVector{3, T}, L::Luneburg)::SVector{3, T}  where {T<:Real}
    if r(x, L) <= L.R
      return (-L.n₀^2 / (n(x, L) * L.R^2)) * (x - L.center)
    else
      return SVector(0.0, 0.0, 0.0)
    end
  end


  """Коэффициент преломления линзы Максвелла"""
  function n(x::SVector{3, T}, L::Maxwell)::Float64 where {T<:Real}
    if r(x, L) <= L.R
        return L.n₀ / (1 + (r(x, L)/L.R)^2)
    else
        return L.n₀
    end
  end

  """Частные производные коэффициента преломления линзы Максвелла"""
  function ∂n(x::SVector{3, T}, L::Maxwell)::SVector{3, T}  where {T<:Real}
    if r(x, L) <= L.R
      return (-2n(x, L)^2 / (L.n₀ * L.R^2)) * (x - L.center)
    else
      return SVector(0.0, 0.0, 0.0)
    end
  end


  """Коэффициент преломления линзы Итона. Линза представляет собой
  диск. Внутренняя окружность радиусом R, а внешняя окружность
  радиусом 2R"""
  function n(x::SVector{3, T}, L::Eaton)::Float64 where {T<:Real}
    if L.R <= r(x, L) <= 2L.R
        return L.n₀ * sqrt( 2L.R / r(x, L) - 1)
    else
        return L.n₀
    end
  end


  """Частные производные коэффициента преломления линзы Итона"""
  function ∂n(x::SVector{3, T}, L::Eaton)::SVector{3, T}  where {T<:Real}
    if L.R <= r(x, L) <= 2L.R
      return -L.n₀ * L.R * (x - L.center) / (sqrt(2L.R / r(x, L) - 1) * r(x, L)^3)
    else
      return SVector(0.0, 0.0, 0.0)
    end
  end


  "Генератор функций для коэффициента преломления"
  function refractive_index_and_its_derivative(L::Lense)
    println("# Установленны коэффициенты преломления для линзы $(L.name)")
    return (
      (x) ->  n(x, L),
      (x) -> ∂n(x, L)
    )
  end

  """Формируем правую часть системы ОДУ уравнения эйконала
  для трехмерного случая в виде канонических уравнений
    x <-> x[1],
    y <-> x[2],
    z <-> x[3],
    p₁<-> x[4],
    p₂<-> x[5],
    p₃<-> x[6].
  """
  function eikonal_equation(L::Lense)

    n, ∂n = refractive_index_and_its_derivative(L)

    function eikonal(v)
      x = SVector{3}(v[1:3])
      p = SVector{3}(v[4:end])
      return SVector(p..., (n(x) * ∂n(x))...)
    end
    println("# Сформировано уравнение эйконала для линзы $(L.name)")
    return (eikonal, n, ∂n)
  end

  """Формируем правую часть системы ОДУ уравнения геодезических кривых
    x <-> x[1],
    y <-> x[2],
    z <-> x[3],
    p₁<-> x[4],
    p₂<-> x[5],
    p₃<-> x[6].
  """
  function geodesic_equation(L::Lense)

    n, ∂n = refractive_index_and_its_derivative(L)

    function geodesic(v)
      x = SVector{3}(v[1:3])
      p = SVector{3}(v[4:end])
      P = (sum(p .* p) + 3 / n(x)^2) / (4n(x))
      return SVector(p..., (P * ∂n(x))...)
    end
    println("# Сформировано уравнение геодезических для линзы $(L.name)")
    return (geodesic, n, ∂n)
  end
end