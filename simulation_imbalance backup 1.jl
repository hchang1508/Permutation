### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 04c7fbc4-4d94-11ed-1276-f58f8d304231
using Random

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
	println(sum(w[1,1:m] .* w[2,1:m]))
	println("In p-value")
	println("Number of treated is: ",m)
	println("Number of treated is: ",n-m)
	
	num_perm = num_perm + 1
	p = 0
	tau_0 = tau_true(w)

	println("samples",samples)
	data_table = zeros(Float64, samples, 5)
	
	if multi_thread==true
		#parallized version
		Threads.@threads for i in 1:samples
			w = w[:, shuffle(1:end)]
          
			#row one is c
			w_data = vcat(w[1, 1:m], w[2,(m+1):n])
			T_w = neyman_estimator(w_data,m,n)
			if abs(T_w - tau_0) >= abs(T_obs - tau_0)
				p += 1
			end
			data_table[i] =  T_w
		end
	else
		#nonparallized version
		for i in 1:samples
			w = w[:, shuffle(1:end)]
			q00_1 = sum((1 .- w[1,1:m]) .* (1 .-w[2,1:m]))
			q00_0 = sum((1 .- w[1,(m+1):n]) .* (1 .-w[2,(m+1):n]))
			q11_1 = sum(w[1,1:m] .* w[2,1:m])
			q11_0 = sum(w[1,(m+1):n] .* w[2,(m+1):n])
			
			
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
	println(data_table[1:2,2:5])
	return p/samples, num_perm,data_table
end

# ╔═╡ fa8efaef-18cf-4f80-adcf-d31da0bd5b01
function create_table(i,j,k,l)
	y  = zeros(Int64, 2, 0)
	y = hcat(y, repeat([1,1],1,i), repeat([1,0],1,j), repeat([0,1],1,k), repeat([0,0],1,l))
	println(size(y))
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
	println("In tau_compatible")
	println("tau0=", t0)
	println("n=", n)
	T_obs = N[1]/m - N[3]/(n-m)

	compatible = false 
	for j in 0:n
		if is_possible(N,j,t0)
			
			v10 = find_v10(N,j,t0)
			v11 = j - v10
			v01 = v10 - t0
			v00 = n - j - v10 + t0
			
			println(v10,v11,v01,v00)
			
			w = create_table(v11,v10,v01,v00)
			p_val, num_perm, data_matrix = p_value(T_obs, w,m,n,samples,num_perm,multi_thread)
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


	max_compatible_tau = N[1] + N[4] 
	min_compatible_tau = -N[2] - N[3]
	
	bin_lower = min_compatible_tau
	bin_upper = max_compatible_tau

	if_compatible = true

	
	counter = 0
	t0=-14
	println("In efficient_perm")
	if_compatible, num_perm =tau_compatible(N,samples,t0,alpha,m,n,num_perm,multi_thread)

	
	return num_perm
	
end

# ╔═╡ 2fca2bd9-8cf3-4cb5-8379-364f81fef446
function perm_test(n11,n10,n01,n00,alpha,m,n,multi_thread,epsilon=0.005,nsample=1000)
	
	#printing basic information
	println("Number of units in the experiment: ",n )
	println("Number of treated with outcome 1: ",n11)
	println("Number of treated with outcome 0: ",n10)
	println("Number of untreated with outcome 1: ",n01)
	println("Number of untreated with outcome 0: ",n00)

	println("Coverage probability of the two-sided interval: ", alpha)

	
	nsample = Int(floor( 1/ (epsilon^2) * log(8 * n * log2(n)/epsilon) ) + 1)
	println("Number of Monte Carlo draws for each permutation test: ", nsample)

	if multi_thread==true
		println("Multi-threading")
	end
	#this variable counts the number of permutation tests conducted
	num_perm ::Int64 = 0

	N_tuple= (n11,n10,n01,n00)

	#decrease alpha by epsilon as equation (35) of the paper
	#alpha = alpha - epsilon
	println("Number of units in treatment: ", m)
	println("Number of units in control: ", n-m)

	num_perm = efficient_perm(N_tuple, nsample, alpha, m,n,0, multi_thread)

	println("")
	#println("95% Confidence Interval:", interval)
	#println("Number of permutation tests: ", num_perm)
end

# ╔═╡ 234558a3-0b15-4963-8b58-b0cc44195f3c
@time perm_test(2,6,8,0,0.05,8,16,false)

# ╔═╡ ceb20b1c-a5d3-426d-b546-e76d3c62972f
@time perm_test(2,6,8,0,0.05,8,8,false)

# ╔═╡ 08a698d2-c4b0-45f9-8a1a-1f8e4f8825bb
@time perm_test(6,4,4,6,0.05,false)

# ╔═╡ 409ef880-13fc-4704-9675-73abdb7f951e
@time perm_test(8,4,5,7,0.05,false)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "fa3e19418881bf344f5796e1504923a7c80ab1ed"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
"""

# ╔═╡ Cell order:
# ╠═d74ce4a0-6f9b-4e63-9728-0208632acf7d
# ╠═04c7fbc4-4d94-11ed-1276-f58f8d304231
# ╠═7f963ab5-958a-4bb1-abf5-82ccfe403b23
# ╠═00b56923-ce85-4a1e-8ed7-68d238d2e07f
# ╠═ea6bc7de-5112-4b5a-9175-01cbcd81ea78
# ╠═234558a3-0b15-4963-8b58-b0cc44195f3c
# ╠═fa8efaef-18cf-4f80-adcf-d31da0bd5b01
# ╠═29577845-c090-4edc-a953-af08d303a592
# ╠═c7432f87-83bf-404d-9618-db05592eb071
# ╠═b4e52381-2725-4cd0-b44c-e17830c5dd52
# ╠═f2cb7f02-6623-4421-b8c0-7e0830d10257
# ╠═8fb19f7e-5d04-4f96-9270-9ffed94147d1
# ╠═966da19d-6713-45de-aa74-4bb743c0b721
# ╠═a2798fc9-24f1-4a76-bb0a-9ae31d09f491
# ╠═2fca2bd9-8cf3-4cb5-8379-364f81fef446
# ╠═04f6c3ce-5204-4641-b555-9dca642a075a
# ╠═ceb20b1c-a5d3-426d-b546-e76d3c62972f
# ╠═08a698d2-c4b0-45f9-8a1a-1f8e4f8825bb
# ╠═409ef880-13fc-4704-9675-73abdb7f951e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
