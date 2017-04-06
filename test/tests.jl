


module AssetTests
    

	using Base.Test, HW_constrained, DataFrames, ForwardDiff
    
    # define "truth"
    truth = DataFrame(a = [0.5, 1.0, 5.0], 
    c = [1.00801, 1.00401, 1.008], 
    omega1 = [-1.41237, -.206197, .758762], 
    omega2 = [.801458, .400729, .0801456], 
    omega3 = [1.60291, .801462, .160291], 
    fval = [-1.20821, -.732819, -.013422])
    
    #define functions
    u(a,c) = -exp(-a*c)
    uprime(a,c) = a*exp(-a*c)

	@testset "testing components" begin

		@testset "finite differences" begin


		end

		@testset "test_finite_diff" begin
		end


		@testset "tests gradient of objective function" begin
            for a in [.5 1.0 1.5]
                d = HW_constrained.data(a)
                x = randn(4)
                obj_t(x) = -u(d["a"],x[1]) - sum(d["pi"]*u(d["a"],sum(x[j+1]*d["z"][s,j] for j in 1:d["n"])) for s in 1:16 ) 
                truegrad = gradient(obj_t, x)
                grad = Vector(4)
                grad[1] = -uprime(d["a"],x[1])
                grad[2] = -sum(d["z"][s,1]*d["pi"]*uprime(d["a"],sum(x[j+1]*d["z"][s,j] for j in 1:d["n"])) for s in 1:16 )
                grad[3] = -sum(d["z"][s,2]*d["pi"]*uprime(d["a"],sum(x[j+1]*d["z"][s,j] for j in 1:d["n"])) for s in 1:16 )
                grad[4] = -sum(d["z"][s,3]*d["pi"]*uprime(d["a"],sum(x[j+1]*d["z"][s,j] for j in 1:d["n"])) for s in 1:16 )

                @test maxabs(truegrad - grad) <1e-3
            end
		end


		@testset "tests gradient of constraint function" begin
          for a in [.5 1 1.5]
                d = HW_constrained.data(a)
                x = randn(4)
                constr_t(x) = -x[1] - sum(d["p"][i]*(x[i+1] - d["e"][i]) for i in 1:d["n"])
                truegrad = gradient(constr_t, x)
                grad = Vector(4)
                grad[1] = -1
                grad[2] = -d["p"][1]
                grad[3] = -d["p"][2]
                grad[4] = -d["p"][3]
                @test maxabs(truegrad - grad) <1e-3
           end
		end
	end

	@testset "testing result of both maximization methods" begin


		@testset "checking result of NLopt maximization" begin
            tno = table_NLopt()
            for i in 1:3, j in 2:6
                @test abs(tno[i,j] - truth[i,j]) < .01
            end
		end


		@testset "checking result of JuMP maximization" begin
            tJu = table_JuMP()
            for i in 1:3, j in 2:6
                @test abs(tJu[i,j] - truth[i,j]) < .01
            end
		end
	end




end



