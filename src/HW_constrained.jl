

# constrained maximization exercises

## portfolio choice problem

module HW_constrained

	using JuMP, NLopt, DataFrames, Ipopt

	export data, table_NLopt, table_JuMP

    function input(prompt::AbstractString="")
       print(prompt)
       return chomp(readline())
    end

    #define u and uprime for NLopt
    u(a,c) = -exp(-a*c)
    uprime(a,c) = a*exp(-a*c)

    # define "truth"
    truth = DataFrame(a = [0.5, 1.0, 5.0], 
    c = [1.00801, 1.00401, 1.008], 
    omega1 = [-1.41237, -.206197, .758762], 
    omega2 = [.801458, .400729, .0801456], 
    omega3 = [1.60291, .801462, .160291], 
    fval = [-1.20821, -.732819, -.013422])

	function data(a=1)
        pi = 1/16
        n = 3
        p = [ 1, 1, 1 ]
        e = [ 2, 0, 0 ]
        z = [repeat([1],inner=16) repeat([.72,.92,1.12,1.32], outer=4) repeat([.86,.96,1.06,1.16], inner=4)]
        return Dict("a" => a, "n" => n, "p" => p, "e" => e, "z" => z, "pi" => pi)
    end
	

	function max_JuMP(a=.5)
        d = data(a)
        m = Model(solver=IpoptSolver())
        @variable(m, c >= 0, start=1.0)
        @variable(m, omega[1:3])

        @NLobjective(m, Max, -exp(-a*c) + sum(-.0625*exp(-a*sum(omega[j]*d["z"][s,j] for j in 1:d["n"])) for s in 1:16 ))
        @NLconstraint(m, c + sum(d["p"][i]*(omega[i] - d["e"][i]) for i in 1:d["n"])==0)

        solve(m)

        cons = getvalue(c)
        maximum = getobjectivevalue(m)
        omegas = getvalue(omega)
        return cons, maximum, omegas
    end

	function table_JuMP()
        tableJuMP = DataFrame(a = Float64[], 
            c = Float64[], omega1 = Float64[], omega2 = Float64[], omega3 = Float64[], fval = Float64[])
        i = 1
        for a in [.5 1 5]
            cons, fval, omegas = max_JuMP(a)
            push!(tableJuMP,[a cons omegas[1] omegas[2] omegas[3] fval])
            i = i+1
        end
        return tableJuMP
    end

	
	function obj(x::Vector,grad::Vector,data::Dict)
        if length(grad) > 0
            grad[1] = -uprime(data["a"],x[1])
            grad[2] = -sum(data["z"][s,1]*data["pi"]*uprime(data["a"],sum(x[j+1]*data["z"][s,j] for j in 1:data["n"])) for s in 1:16 )
            grad[3] = -sum(data["z"][s,2]*data["pi"]*uprime(data["a"],sum(x[j+1]*data["z"][s,j] for j in 1:data["n"])) for s in 1:16 )
            grad[4] = -sum(data["z"][s,3]*data["pi"]*uprime(data["a"],sum(x[j+1]*data["z"][s,j] for j in 1:data["n"])) for s in 1:16 )
        end
        return -u(data["a"],x[1]) - sum(data["pi"]*u(data["a"],sum(x[j+1]*data["z"][s,j] for j in 1:data["n"])) for s in 1:16 )
    end

	function constr(x::Vector,grad::Vector,data::Dict)
        if length(grad) > 0 
            grad[1] = -1
            grad[2] = -data["p"][1]
            grad[3] = -data["p"][2]
            grad[4] = -data["p"][3]
        end
        return -x[1] - sum(data["p"][i]*(x[i+1] - data["e"][i]) for i in 1:data["n"])
    end

    function constr2(x::Vector,grad::Vector,data::Dict)
        if length(grad) > 0 
            grad[1] = 1
            grad[2] = data["p"][1]
            grad[3] = data["p"][2]
            grad[4] = data["p"][3]
        end
        return x[1] + sum(data["p"][i]*(x[i+1] - data["e"][i]) for i in 1:data["n"])
    end

	function max_NLopt(a=0.5)
		opt = Opt(:LD_MMA, 4)
        # set bounds and tolerance
        lower_bounds!(opt, [0, -Inf, -Inf, -Inf])
        #upper_bounds!(opt, [2, 2, 2, 2])
        xtol_rel!(opt,1e-4)
        d = data(a)
        obj2(x,g) = obj(x,g,d)
        # define objective function
        min_objective!(opt, obj2)
        # define constraints
        # notice the anonymous function
        inequality_constraint!(opt, (x,g) -> constr(x,g,d), 1e-8)
        inequality_constraint!(opt, (x,g) -> constr2(x,g,d), 1e-8)
        # call optimize
        return optimize(opt, [1.0, 1.0, 0.0, 0.0])
	end

	function table_NLopt()
        tableNLopt = DataFrame(a = Float64[], 
            c = Float64[], omega1 = Float64[], omega2 = Float64[], omega3 = Float64[], fval = Float64[])
        i = 1
        for a in [.5 1 5]
            #d = data(a)
            fval, xvals, res = max_NLopt(a)
            push!(tableNLopt,[a xvals[1] xvals[2] xvals[3] xvals[4] -fval])
            i = i+1
        end
        return tableNLopt
    end

	# function `f` is for the NLopt interface, i.e.
	# it has 2 arguments `x` and `grad`, where `grad` is
	# modified in place
	# if you want to call `f` with more than those 2 args, you need to
	# specify an anonymous function as in
	# other_arg = 3.3
	# test_finite_diff((x,g)->f(x,g,other_arg), x )
	# this function cycles through all dimensions of `f` and applies
	# the finite differencing to each. it prints some nice output.
	function test_finite_diff(f::Function,x::Vector{Float64},tol=1e-6)
	
	end

	# do this for each dimension of x
	# low-level function doing the actual finite difference
	function finite_diff(f::Function,x::Vector)
		
	end

	function runAll()
		println("running tests:")
		include("test/runtests.jl")
		println("")
		println("JumP:")
		display(table_JuMP())
		println("")
		println("NLopt:")
		display(table_NLopt())
        println("Truth:")                                                                                                                                
        display(truth)                                                                                                              
		ok = input("enter y to close this session.")
		if ok == "y"
			quit()
		end
	end


end


