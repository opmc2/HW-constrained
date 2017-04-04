


module AssetTests

	using Base.Test, HW_constrained, DataFrames

    truth = DataFrame(a = [0.5, 1.0, 5.0], 
    c = [1.00801, 1.00401, 1.008], 
    omega1 = [-1.41237, -.206197, .758762], 
    omega2 = [.801458, .400729, .0801456], 
    omega3 = [1.60291, .801462, .160291], 
    fval = [-1.20821, -.732819, -.013422])

	@testset "testing components" begin

		@testset "finite differences" begin


		end

		@testset "test_finite_diff" begin
		end


		@testset "tests gradient of objective function" begin
		end


		@testset "tests gradient of constraint function" begin
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



