using Test
using PoissonBinomialTest
p = Probability
@test p(0.5, 1) == p(0.5, 0)
@test p(0.3) < p(0.9)
@test p(0.3) < p(0.4)
@test p(0.8) < p(0.9)
@test ! (p(0.9) < p(0.3))
dist = PBDist([0.3, 0.9])
edist = expand(dist)
@test isapprox(edist, [0.07, 0.66, 0.27])
@test isapprox(collect(float.(pbtest(dist, 0))), [0.00, 0.93])
@test isapprox(collect(float.(pbtest(dist, 1))), [0.07, 0.27])
@test isapprox(collect(float.(pbtest(dist, 2))), [0.73, 0.00])
