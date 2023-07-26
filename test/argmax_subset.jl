using BenchmarkTools
using BenchmarkPlots
using StatsPlots
using Test
using Random: bitrand

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 20
N = 10_000

preservingargmax = Dict{Symbol,Function}()

preservingargmax[:Inf] = (x, b) -> argmax(x .- Inf .* .!(b))
preservingargmax[:findall] = (x, b) -> findall(b)[argmax(x[b])]
preservingargmax[:findall_views] = (x, b) -> first(@view findall(b)[argmax(@view x[b])])
preservingargmax[:findnext] = (x, b) -> begin
    c = argmax(x[b])
    i0 = findfirst(identity, b)
    for _ in 1:c-1
        i0 = findnext(identity, b, i0 + 1)
    end
    return i0
end
preservingargmax[:findnext_views] = (x, b) -> begin
    c = argmax(view(x, b))
    i0 = findfirst(identity, b)
    for _ in 1:c-1
        i0 = findnext(identity, b, i0 + 1)
    end
    return i0
end
preservingargmax[:eachindex] = (x, b) -> eachindex(x)[b][argmax(x[b])]
preservingargmax[:eachindex_views] = (x, b) -> first(view(view(eachindex(x),b),argmax(view(x,b))))
preservingargmax[:sortperm] = (x, b) -> argmax(sortperm(x) .* b)
preservingargmax[:argmax] = (x, b) -> argmax(i -> b[i] ? x[i] : -Inf, eachindex(x))
preservingargmax[:enumerate] = (x, b) -> first(argmax(t -> t[2][1] ? t[2][2] : -Inf, enumerate(zip(b, x))))
preservingargmax[:zip] = (x, b) -> first(argmax(t -> t[2] ? t[3] : -Inf, zip(eachindex(x), b, x)))
preservingargmax[:loop] = (x, b) -> begin
    maxx = -Inf
    maxi = 0
    for (xi, bi, i) in zip(x, b, eachindex(x))
        (bi && xi > maxx) || continue
        maxx = xi
        maxi = i
    end
    return maxi
end
preservingargmax[:branchless] = (x, b) -> begin
    maxx = -Inf
    maxi = 0
    for (xi, bi, i) in zip(x, b, eachindex(x))
        maxi = (bi && xi > maxx)*i + !(bi && xi > maxx)*maxi
        maxx = (bi && xi > maxx)*xi + !(bi && xi > maxx)*maxx
    end
    return maxi
end
preservingargmax[:to_index] = (x, b) -> (i=Base.to_index(b); Vector{Int}(i)[argmax(x[i])])
preservingargmax[:typemin] = (x, b) -> (y = copy(x); y[.!b] .= typemin(eltype(x)); argmax(y))
preservingargmax[:parentindices] = (x, b) -> (v = view(x, b); first(parentindices(v))[argmax(v)])

labels = [
        :Inf,
        :findall,
        :findall_views,
        :findnext,
        :findnext_views,
        :eachindex,
        :eachindex_views,
        :sortperm,
        :argmax,
        :enumerate,
        :zip,
        :loop,
        :branchless,
        :to_index,
        :typemin,
        :parentindices,
    ]

bg = BenchmarkGroup()

x_1 = [0, 2, 1, -3, 0]
b_1 = [true, false, true, true, false]
x_2 = collect(1:100)
b_2 = falses(100); b_2[[1, 3, 10, 100]] .= true
@testset "correctness, $k" for (k, f) in preservingargmax
    @test f(x_1, b_1) == 3
    @test f(x_2, b_2) == 100
    bg[k] = @benchmarkable $f(x, b) setup = (x = randn($N); b = bitrand($N))
end

res = run(bg)


plot(
    res,
    labels,
    ;
    yscale=:log10,
    xrotation=30,
)

savefig("benchmark.png")
