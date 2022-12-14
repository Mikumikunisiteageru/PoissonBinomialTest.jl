module PoissonBinomialTest

export PBDist, Probability, expand, pbtest

export complement, log, logit

const ALMOST = 0.999

struct Probability{T <: AbstractFloat}
	off::T
	ref::Bool
	function Probability(off::T, ref::Bool=false) where {T <: AbstractFloat}
		0 <= off <= 1 || error("Range error!")
		off > ALMOST && @warn "Low precision when transforming $(off)!"
		off == 0.5 && return new{T}(off, false)
		return off > 0.5 ? new{T}(1-off, !ref) : new{T}(off, ref)
	end
end
Probability(off::T, ref::Integer) where {T <: AbstractFloat} =
	Probability(off, Bool(ref)) # Parameter `ref` has to be either zero or one

function Base.show(io::IO, prob::Probability)
	print(io, prob.ref ? "1 - $(repr(prob.off))" : "0 + $(repr(prob.off))")
end

function Base.isless(a::Probability, b::Probability)
	iszero(a.ref) && iszero(b.ref) && return a.off < b.off
	isone(a.ref)  && isone(b.ref)  && return b.off < a.off
	return (iszero(a.ref) ? (<) : (>))(a.off + b.off, 1)
end

Base.float(p::Probability) = p.ref ? 1 - p.off : p.off
Base.isapprox(a::Probability, b::Probability) = 
	isapprox(1-(1-(float(a))), 1-(1-(float(b))))
Base.length(p::Probability) = 1
Base.iterate(p::Probability) = (p, nothing)
Base.iterate(p::Probability, nothing) = nothing

# Base.:+(p::Probability, x::Number) = p.ref ? x + 1 - p.off : x + p.off
# Base.:+(x::Number, p::Probability) = p.ref ? x + 1 - p.off : x + p.off
# Base.:-(p::Probability, x::Number) = p.ref ? 1 - x - p.off : p.off - x
# Base.:-(x::Number, p::Probability) = p.ref ? x - 1 + p.off : x - p.off

complement(p::Probability) = Probability(p.off, !p.ref)
Base.log(p::Probability) = p.ref ? log1p(-p.off) : log(p.off)
logit(p::AbstractFloat) = log(p) - log(1 - p)
logit(p::Probability) = log(p) - log(complement(p))

struct PBDist{T <: AbstractFloat}
	n::Int
	pp::Vector{T}
	PBDist(pp::Vector{T}) where {T <: AbstractFloat} =
		new{T}(length(pp), pp)
end

function expand(dist::PBDist{T}) where {T <: AbstractFloat}
	slots = zeros(T, dist.n + 1)
	slots[1] = 1
	cnt = 1
	for x = dist.pp
		y = 1 - x
		iszero(x) && continue
		cnt += 1
		for i = cnt:-1:2
			slots[i] = slots[i] * y + slots[i-1] * x
		end
		slots[1] *= y
	end
	slots
end

function pbtest(dist::PBDist{<:AbstractFloat}, x::Int)
	@assert 0 <= x <= dist.n
	slots = expand(dist)
	left  = sum(slots[1:x])
	mid   = slots[x+1]
	right = sum(reverse(slots[x+2:end]))
	p_left  = left  < 0.5 ? Probability(left) : Probability(mid + right, true)
	p_right = right < 0.5 ? Probability(right) : Probability(mid + left, true)
	p_left, p_right
end
pbtest(pp::Vector{<:AbstractFloat}, x::Int) = pbtest(PBDist(pp), x)

function quantile(dist::PBDist{<:AbstractFloat}, x::Int)
	@assert 0 <= x <= dist.n
	slots = expand(dist)
	mleft  = sum(slots[1:x]) + slots[x+1] / 2
	mright = sum(reverse(slots[x+2:end])) + slots[x+1] / 2
	mleft < 0.5 ? Probability(mleft) : Probability(mright, true)
end
quantile(pp::Vector{<:AbstractFloat}, x::Int) = quantile(PBDist(pp), x)

end # module
