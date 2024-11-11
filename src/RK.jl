module RK

using Base.Iterators: countfrom, takewhile
using StaticArrays: SVector
"""
Метод Рунге--Кутты с постоянным шагом интегрирования. Результаты вычислений и
моменты времени записываются в массивы и возвращаются.

```RKp6n1(func, x_0, h, start, stop) -> (T, X)```

# Аргументы

* `func::Function`: правая часть системы ОДУ.
* `x_0::SVector`: начальное значение `x(t)`.
* `h::Float64`: шаг метода.
* `start::Float64`: начальная точка временного интервала.
* `stop::Float64`: конечная точка временного интервала.

# Возвращаемые значения

* `T::Array{Float64, 1}`: значения времени из интервала `[start, stop]`.
* `X::Array{Float64, EQN}`: численное решение в виде массива из `N` строк и  `EQN` столбцов.
"""
function RKp6n1(func::Function, x_0::SVector, h::Float64, start::Float64, stop::Float64)

  # Массив, содержащий точки временной сетки
  T = collect(takewhile(<=(stop), countfrom(start, h)))

  EQN = length(x_0)
  N = length(T)

  # В массиве X N строк и EQN столбцов
  X = Matrix{Float64}(undef, N, EQN)

  x = SVector{EQN}(x_0)

  for i ∈ 1:N
    X[i,:] = x
    k1 = func(x)
    k2 = func(x + h*(1//2*k1))
    k3 = func(x + h*(2//9*k1 + 4//9*k2))
    k4 = func(x + h*(7//36*k1 + 2//9*k2 + -1//12*k3))
    k5 = func(x + h*(-35//144*k1 + -55//36*k2 + 35//48*k3 + 15//8*k4))
    k6 = func(x + h*(-1//360*k1 + -11//36*k2 + -1//8*k3 + 1//2*k4 + 1//10*k5))
    k7 = func(x + h*(-41//260*k1 + 22//13*k2 + 43//156*k3 + -118//39*k4 + 32//195*k5 + 80//39*k6))

    x = x + h*(13//200*k1 + 11//40*k3 + 11//40*k4 + 4//25*k5 + 4//25*k6 + 13//200*k7)

  end
  return (T, X)
end

"""Модифицированный вариант RKp6n1, который возвращает не кортеж, а массив
```RKp6n1(func, x_0, h, start, stop) -> res```
"""
function RKp6n1_v2(; func::Function, x_0::SVector, h::Float64, start::Float64, stop::Float64)::Matrix{Float64}

  # Массив, содержащий точки временной сетки
  T = collect(takewhile(<=(stop), countfrom(start, h)))

  EQN = length(x_0)
  N = length(T)
  
  # В массиве X N строк и EQN+1 столбцов
  # где в первом столбце записывается время, а в остальных координаты
  # соответствующие этому времени
  res = Matrix{Float64}(undef, N, EQN+1)
  res[:, 1] = T
  x = SVector{EQN}(x_0)

  for i ∈ 1:N
    res[i,2:end] = x
    k1 = func(x)
    k2 = func(x + h*(1//2*k1))
    k3 = func(x + h*(2//9*k1 + 4//9*k2))
    k4 = func(x + h*(7//36*k1 + 2//9*k2 + -1//12*k3))
    k5 = func(x + h*(-35//144*k1 + -55//36*k2 + 35//48*k3 + 15//8*k4))
    k6 = func(x + h*(-1//360*k1 + -11//36*k2 + -1//8*k3 + 1//2*k4 + 1//10*k5))
    k7 = func(x + h*(-41//260*k1 + 22//13*k2 + 43//156*k3 + -118//39*k4 + 32//195*k5 + 80//39*k6))

    x = x + h*(13//200*k1 + 11//40*k3 + 11//40*k4 + 4//25*k5 + 4//25*k6 + 13//200*k7)

  end
  return res
end

end #module RK
