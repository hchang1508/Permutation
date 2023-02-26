### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 04c7fbc4-4d94-11ed-1276-f58f8d304231
using Random

# ╔═╡ cd732c96-81f5-4356-8736-5685f5bad59d
using Distributions

# ╔═╡ d74ce4a0-6f9b-4e63-9728-0208632acf7d
#This file creates simulation results in Table 1, Column (Algorithm-4.3) in the paper.

# ╔═╡ 7f963ab5-958a-4bb1-abf5-82ccfe403b23
function neyman_estimator(data,m,n)
	T = sum(data[1:m])/m - sum(data[m+1:n])/(n-m)
end

# ╔═╡ 00b56923-ce85-4a1e-8ed7-68d238d2e07f
function tau_true(table)
	return (sum(table[1,:]) - sum(table[2,:]))/length(table[1,:])
end

# ╔═╡ ea6bc7de-5112-4b5a-9175-01cbcd81ea78
function p_value(T_obs, w, m,n,samples,num_perm,multi_thread)
	
	#println(num_perm)
	#println(w)
	#open("file.txt","a") do io
   	#	println(io,"permutation test performed")
	#end

	num_perm = num_perm + 1
	p = 0
	tau_0 = tau_true(w)

	data_table = zeros(Float64, samples, 5)
	
	if multi_thread==true
		#parallized version
		Threads.@threads for i in 1:samples
			
			w = w[:, shuffle(1:end)]

			#number of 00 and 11 assigned to treatment and control
			q00_1 = sum((1 .- w[1,1:m]) .* (1 .-w[2,1:m]))
			q00_0 = sum((1 .- w[1,(m+1):n]) .* (1 .-w[2,(m+1):n]))
			q11_1 = sum(w[1,1:m] .* w[2,1:m])
			q11_0 = sum(w[1,(m+1):n] .* w[2,(m+1):n])

			println("q00_1: ", q00_1)
			println("q00_0: ", q00_0)
			println("q11_1: ", q11_1)
			println("q11_0: ", q11_0)
			
			#realized data
			w_data = vcat(w[1, 1:m], w[2,(m+1):n])
			T_w = neyman_estimator(w_data,m,n)

			
			if abs(T_w - tau_0) >= abs(T_obs - tau_0)
				p += 1
			end
			data_table[i,1] =  T_w
			data_table[i,2] =  q00_1
			data_table[i,3] =  q00_0
			data_table[i,4] =  q11_1
			data_table[i,5] =  q11_0

		end		
	else
		#nonparallized version
		for i in 1:samples
			
			w = w[:, shuffle(1:end)]

			#number of 00 and 11 assigned to treatment and control			
			q00_1 = sum((1 .- w[1,1:m]) .* (1 .-w[2,1:m]))
			q00_0 = sum((1 .- w[1,(m+1):n]) .* (1 .-w[2,(m+1):n]))
			q11_1 = sum(w[1,1:m] .* w[2,1:m])
			q11_0 = sum(w[1,(m+1):n] .* w[2,(m+1):n])
			
			#w_data			
			w_data = vcat(w[1, 1:m], w[2,(m+1):n])
			T_w = neyman_estimator(w_data,m,n)

			
			if abs(T_w - tau_0) >= abs(T_obs - tau_0)
				p += 1
			end
			data_table[i,1] =  T_w
			data_table[i,2] =  q00_1
			data_table[i,3] =  q00_0
			data_table[i,4] =  q11_1
			data_table[i,5] =  q11_0

		end		
	end
	return p/samples, num_perm, data_table
end

# ╔═╡ fa8efaef-18cf-4f80-adcf-d31da0bd5b01
function create_table(i,j,k,l)
	y  = zeros(Int64, 2, 0)
	y = hcat(y, repeat([1,1],1,i), repeat([1,0],1,j), repeat([0,1],1,k), repeat([0,0],1,l))
	return y
end

# ╔═╡ 29577845-c090-4edc-a953-af08d303a592
function check_possible(N,V)
	return ( max(0, N[1] - V[2], V[1] - N[3], V[1] + V[3] - N[2] - N[3]) 
		<= min(V[1], N[1], V[1] + V[3] - N[3], sum(N) - V[2] - N[3] - N[2]) )
end

# ╔═╡ c7432f87-83bf-404d-9618-db05592eb071
function median(a::Int64, b::Int64)
	return floor(Int64,(a+b)/2)
end

# ╔═╡ b4e52381-2725-4cd0-b44c-e17830c5dd52
function is_possible(N,j,t0)

	#This function checks whether the potential outcome table, indexed by j, is compatible with actual data, N, and the hypothesis to be tested, t0
	
	flag1 = (j >= t0 + N[3]) && (j >= N[1]) && (sum(N)>=j+N[2]) && (N[1] + t0 + N[2] + N[3] >= j)
	flag2 = max(t0, j - N[1] - N[3], N[1] + N[3] + t0 -j) <= min(j, N[1] + N[4], N[2] + N[3] + t0, sum(N) + t0 -j)

	return flag1 && flag2
end

# ╔═╡ f2cb7f02-6623-4421-b8c0-7e0830d10257
function find_v10(N,j,t0)
	return max(t0, j - N[1] - N[3], N[1] + N[3] + t0 -j,0)
end 

# ╔═╡ 8fb19f7e-5d04-4f96-9270-9ffed94147d1
function find_v10_upper(N,j,t0)
	return min(j, N[1] + N[4], N[2] + N[3] + t0, sum(N) + t0 - j)
end

# ╔═╡ 77e0305b-84c6-47a5-a736-8ce2b61f12c9
function update_pval(data_matrix,t0,samples,T_obs,n)

	println("update_pval")
	p =  abs.(data_matrix[:,1] .- t0/n) .> abs(T_obs-t0/n)
 	p = Int.(p)
	p_sum = sum(p)/samples
	
	return p_sum


end

# ╔═╡ 968ca6e1-651d-4782-a275-d89bb226593f
function update(data_matrix,samples,m,n)

	q00 = rand(Uniform(0,1),samples)
	q11 = rand(Uniform(0,1),samples)

	
	p11=data_matrix[:,3]./(data_matrix[:,2] .+ data_matrix[:,3] ) 
	p00=data_matrix[:,5]./(data_matrix[:,4] .+ data_matrix[:,5] ) 

	draw00 = Int.(q00 .< p00)
	draw11 = Int.(q11 .< p11)

	update_HT = draw00 .* -1/(n-m) .+ draw11 .* 1/(n-m)
	update = hcat(update_HT, (-1).^ (1 .-draw00), (-1).^ -draw00,(-1).^ (1 .-draw11),(-1).^ draw11)

	data_matrix = data_matrix .+ update
	println("mean is: ", mean(data_matrix[:,1]))
	return data_matrix
end

# ╔═╡ 966da19d-6713-45de-aa74-4bb743c0b721
function tau_compatible(N,samples,t0,alpha,m,n,num_perm,multi_thread) #here t0 is an integer
	
	###inputs
	# T_obs: estimator using the actual data
	# samples: number of Monte Carlo draws
	# t0: hypothesized value to be tested
	# m: number of treated
	# n: number of experimental units
	# alpha: coverage of the interval
	# num_perm: (global variable) to keep track of the number of permutation tests
	# multi_thread: whether the MC procedures are paralleled or not
	###

	#estimated ate
	T_obs = N[1]/m - N[3]/(n-m)

	compatible = false 

	#second line search (first one is tau0)
	for j in 0:n
		
		if is_possible(N,j,t0)

			println("Second loop: ", j)
            #find lowest possible v10
			v10 = find_v10(N,j,t0)
			v11 = j - v10
			v01 = v10 - t0
			v00 = n - j - v10 + t0

			#create potential outcome table
			w = create_table(v11,v10,v01,v00)

			#endpoint permutation test -> save data
			p_val, num_perm, data_matrix = p_value(T_obs, w,m,n,samples,num_perm,multi_thread)

			#second dimension
			
			for k in 1:n

				v10 = v10 + 1 
				v11 = v11 - 1
				v01 = v01 + 1 
		    	v00 = v00 - 1

				V= [v11,v10,v01,v00]
				#print("Data Table")
				#println(V)
				
				if check_possible(N,V)
					println("Third loop: ", k)
					println("It is a possible data table.")
					println(V)					
					data_matrix= update(data_matrix,samples,m,n)

					#Calculate p_val based oon the new results
					p_val_new = update_pval(data_matrix,t0,samples,T_obs,n)

                    #If larger, update
					println("Old p_val:, ",p_val)
					println("New p_val:, ",p_val_new)

					p_val = max(p_val,p_val_new)
				end
			end
			println("Final p value is ",p_val)
			println(alpha)
			if p_val >= alpha
				compatible = true
				break 
			end
		end
	end

	return compatible, num_perm
end

# ╔═╡ a2798fc9-24f1-4a76-bb0a-9ae31d09f491
function efficient_perm(N, samples, alpha,m,n,num_perm,multi_thread)


	
	max_compatible_tau = N[1] + N[4] #maximum number of people with effects when treated: 1 when treated and 0 if in control
	min_compatible_tau = -N[2] - N[3] #minimum number of people with effects when treated: 0 when treated and 1 if in control
	
	bin_lower = min_compatible_tau
	bin_upper = max_compatible_tau

	if_compatible = true

	comptatible_effects=Array{Int64}(undef,0)
	counter = 0
	
	for t0 in min_compatible_tau:max_compatible_tau
		println("----------------------------")
		println("tau0=", t0)
		if_compatible, num_perm =tau_compatible(N,samples,t0,alpha,m,n,num_perm,multi_thread)

		println("tau0 is compatible?", if_compatible)
		println("----------------------------")
		
		if if_compatible==true
			comptatible_effects=push!(comptatible_effects,t0)
		end
	end
	println(num_perm)
	
	println(comptatible_effects)
	
	return num_perm,comptatible_effects
	
end

# ╔═╡ 2fca2bd9-8cf3-4cb5-8379-364f81fef446
function perm_test(n11,n10,n01,n00,alpha,m,n,multi_thread,epsilon=0.001,nsample=1000)
	
	#printing basic information
	println("Number of units in the experiment: ",n )
	println("Number of treated with outcome 1: ",n11)
	println("Number of treated with outcome 0: ",n10)
	println("Number of untreated with outcome 1: ",n01)
	println("Number of untreated with outcome 0: ",n00)

	println("Coverage probability of the two-sided interval: ", alpha)

	
	#nsample = Int(floor( 1/ (epsilon^2) * log(4 * n ^3/epsilon) ) + 1)

	nsample =100000
	println("Number of Monte Carlo draws for each permutation test: ", nsample)

	if multi_thread==true
		println("Multi-threading")
	end
	#this variable counts the number of permutation tests conducted
	num_perm ::Int64 = 0

	N_tuple= (n11,n10,n01,n00)
	println(N_tuple)
	#decrease alpha by epsilon as equation (35) of the paper
	alpha = alpha - epsilon
	println("Number of units in treatment: ", m)
	println("Number of units in control: ", n-m)

	num_perm,comptatible_effects = efficient_perm(N_tuple, nsample, alpha, m,n,0, multi_thread)

	#println("")
	#println("95% Confidence Interval:", interval)
	println("Number of permutation tests: ", num_perm)
end

# ╔═╡ e94f0357-57c3-4afa-a601-74c7bd21def8
perm_test(2, 6, 8, 0,0.05,8,16,true)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
Distributions = "~0.25.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "1452cf99a06750e548f6e193ccacba3c65e7087d"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9a782b47da6ee4cb3d041764e0c6830469105984"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.83"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "d3ba08ab64bdfd27234d3f61956c966266757fe6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.7"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "637b58b3c037d3877f263418de820920b47ceeb5"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "786efa36b7eff813723c4849c90456609cf06661"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "53cd758b96903d556e96f11b8cd2169c7e9f08af"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.2.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═d74ce4a0-6f9b-4e63-9728-0208632acf7d
# ╠═04c7fbc4-4d94-11ed-1276-f58f8d304231
# ╠═cd732c96-81f5-4356-8736-5685f5bad59d
# ╠═7f963ab5-958a-4bb1-abf5-82ccfe403b23
# ╠═00b56923-ce85-4a1e-8ed7-68d238d2e07f
# ╠═ea6bc7de-5112-4b5a-9175-01cbcd81ea78
# ╠═fa8efaef-18cf-4f80-adcf-d31da0bd5b01
# ╠═29577845-c090-4edc-a953-af08d303a592
# ╠═c7432f87-83bf-404d-9618-db05592eb071
# ╠═b4e52381-2725-4cd0-b44c-e17830c5dd52
# ╠═f2cb7f02-6623-4421-b8c0-7e0830d10257
# ╠═8fb19f7e-5d04-4f96-9270-9ffed94147d1
# ╠═966da19d-6713-45de-aa74-4bb743c0b721
# ╠═77e0305b-84c6-47a5-a736-8ce2b61f12c9
# ╠═968ca6e1-651d-4782-a275-d89bb226593f
# ╠═a2798fc9-24f1-4a76-bb0a-9ae31d09f491
# ╠═b0734b78-590a-4c12-a817-1369a6549056
# ╠═e90d8401-2933-4f18-be36-476d8151e29f
# ╠═e94f0357-57c3-4afa-a601-74c7bd21def8
# ╠═7a950e53-c0d4-4ed0-bfae-b0e837161687
# ╠═ea588022-0c1f-4f87-a5e3-c0325f652f42
# ╠═c7591e11-3a03-440c-86e3-d40807be2934
# ╠═2fca2bd9-8cf3-4cb5-8379-364f81fef446
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
