### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 04c7fbc4-4d94-11ed-1276-f58f8d304231
using Random

# ╔═╡ 2ccb99d4-8be1-4db0-8d93-a524c7514278
const alpha = .1

# ╔═╡ 7f963ab5-958a-4bb1-abf5-82ccfe403b23
function neyman_estimator(data)
	m = floor(Int64,length(data)/2)
	T = sum(data[1:m])/m - sum(data[m+1:2*m])/m
end

# ╔═╡ 00b56923-ce85-4a1e-8ed7-68d238d2e07f
function tau_true(table)
	return (sum(table[1,:]) - sum(table[2,:]))/length(table[1,:])
end

# ╔═╡ ea6bc7de-5112-4b5a-9175-01cbcd81ea78
function p_value(T_obs, w, samples)
	open("file.txt","a") do io
   		println(io,"permutation test performed")
	end
	p = 0
	m = floor(Int64, length(w[1,:])/2)
	tau_0 = tau_true(w)
	for i in 1:samples
		w = w[:, shuffle(1:end)]
		w_data = vcat(w[1, 1:m], w[2,m+1:2*m])
		T_w = neyman_estimator(w_data)
		if abs(T_w - tau_0) >= abs(T_obs - tau_0)
			p += 1
		end
	end
	return p/samples
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

# ╔═╡ 923d9586-6453-4d6d-8438-e352c118cda7
function perm(N,samples) #N is [n11,n10,n01,n00] in the notation of the note
	lower = 2.0 # arbitrary values to initialize
	upper = -2.0
	T_obs_rescaled = 2*(N[1] - N[3]) #rescaled by n compared to paper
	max_compatible_tau = N[1] + N[4] #rescaled by n compared to paper
	min_compatible_tau = -N[2] - N[3] 
	n = sum(N)
	#@show T_obs_rescaled/n

	#first find the upper bound
	for x in reverse(min_compatible_tau:max_compatible_tau)
		for i in 0:n, j in 0:n #i is v11, j is v10
			k = j - x #v01
			if k < 0 continue end
			l = n - i - j - k #v00
			if l < 0 continue end
			if check_possible(N, [i,j,k,l]) == false continue end
			w = create_table(i,j,k,l)
			if p_value(T_obs_rescaled/n, w, samples) >= alpha
				upper = x/n
				#@show (i,j,k,l)
				break
			end
		end
		if upper > -1
			break
		end
	end

	for x in min_compatible_tau:max_compatible_tau
		for i in 0:n, j in 0:n #i is v11, j is v10
			k = j - x #v01
			if k < 0 continue end
			l = n - i - j - k #v00
			if l < 0 continue end
			if check_possible(N, [i,j,k,l]) == false continue end
			w = create_table(i,j,k,l)
			#println(i,j,k,l)
			if p_value(T_obs_rescaled/n, w, samples) >= alpha
				lower = x/n
				#@show (i,j,k,l)
				break
			end
		end
		if lower < 1
			break
		end
	end

	return n*[lower,upper]
end
	

# ╔═╡ c7432f87-83bf-404d-9618-db05592eb071
function median(a::Int64, b::Int64)
	return floor(Int64,(a+b)/2)
end

# ╔═╡ b4e52381-2725-4cd0-b44c-e17830c5dd52
function is_possible(N,j,t0)
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
function tau_compatible(N,samples,t0) #here t0 is an integer
	compatible = false 
	T_obs_rescaled = 2*(N[1] - N[3])
	n = sum(N)
	for j in 0:n
		if is_possible(N,j,t0)
			v10 = find_v10(N,j,t0)
			v11 = j - v10
			v01 = v10 - t0
			v00 = n - j - v10 + t0
			w = create_table(v11,v10,v01,v00)
			if p_value(T_obs_rescaled/n, w, samples) >= alpha
				compatible = true
				break 
			end

			if v10 == 0 && v01 == 0
				v10 = 1
				v11 = j - v10
				v01 = v10 - t0
				v00 = n - j - v10 + t0
				if v11 >= 0 && v11 >= 0 && v00 >= 0 && v10 <= find_v10_upper(N,j,t0)
					if p_value(T_obs_rescaled/n, w, samples) >= alpha
					compatible = true
					break 
					end
				end
			end
		end
	end

	return compatible
end

# ╔═╡ a2798fc9-24f1-4a76-bb0a-9ae31d09f491
function efficient_perm(N,samples)
	T_obs_rescaled = 2*(N[1] - N[3])
	#m = floor(Int64, sum(N)/2)
	#T_obs = T_obs_rescaled/m
	max_compatible_tau = N[1] + N[4] #rescaled by n
	min_compatible_tau = -N[2] - N[3]
	#upper = 1
	#lower = -1

	bin_lower = T_obs_rescaled
	bin_upper = max_compatible_tau
	#println(bin_lower,bin_upper)

	counter = 0
	while counter < 50
		#println(counter)
		counter += 1 
		if bin_lower == bin_upper 
			break
		end

		if bin_lower == bin_upper - 1
			if tau_compatible(N,samples,bin_upper) == true
				bin_lower = bin_upper
				break
			else 
				bin_upper = bin_lower
				break
			end
		end
			

		t0 = median(bin_lower, bin_upper)
		
		compatible = tau_compatible(N,samples,t0)
		
		if compatible == true
			bin_lower = t0
		end

		if compatible == false
			bin_upper = t0
		end
		
	end

	n = sum(N)
	upper = bin_lower/n

	bin_lower = min_compatible_tau
	bin_upper = T_obs_rescaled

	counter = 0
	while counter < 50
		#println(counter)
		counter += 1 
		if bin_lower == bin_upper 
			break
		end

		if bin_lower == bin_upper - 1
			if tau_compatible(N,samples,bin_lower) == true
				bin_upper = bin_lower
				break
			else 
				bin_lower = bin_upper
				break
			end
		end
			

		t0 = median(bin_lower, bin_upper)
		
		compatible = tau_compatible(N,samples,t0)
		
		if compatible == true
			bin_upper = t0
		end

		if compatible == false
			bin_lower = t0
		end
		
	end

	lower = bin_lower/n
	return n*[lower,upper]
	
end

# ╔═╡ 70a93020-ddc9-4b9c-b865-0084f26962fe
efficient_perm((40,160,180,20), 1000)

# ╔═╡ 23e787d0-c6c8-4b12-8fff-b01caff74a2e
perm((4,4,4,4),100)

# ╔═╡ 10f4a390-d6aa-4de9-a882-aa68e48064d9
html"""
<style>
  pluto-helpbox {
    display: none;
  }
</style>
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.1"
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
# ╠═04c7fbc4-4d94-11ed-1276-f58f8d304231
# ╠═2ccb99d4-8be1-4db0-8d93-a524c7514278
# ╠═7f963ab5-958a-4bb1-abf5-82ccfe403b23
# ╠═00b56923-ce85-4a1e-8ed7-68d238d2e07f
# ╠═ea6bc7de-5112-4b5a-9175-01cbcd81ea78
# ╠═fa8efaef-18cf-4f80-adcf-d31da0bd5b01
# ╠═29577845-c090-4edc-a953-af08d303a592
# ╠═923d9586-6453-4d6d-8438-e352c118cda7
# ╠═c7432f87-83bf-404d-9618-db05592eb071
# ╠═b4e52381-2725-4cd0-b44c-e17830c5dd52
# ╠═f2cb7f02-6623-4421-b8c0-7e0830d10257
# ╠═8fb19f7e-5d04-4f96-9270-9ffed94147d1
# ╠═966da19d-6713-45de-aa74-4bb743c0b721
# ╠═a2798fc9-24f1-4a76-bb0a-9ae31d09f491
# ╠═70a93020-ddc9-4b9c-b865-0084f26962fe
# ╠═23e787d0-c6c8-4b12-8fff-b01caff74a2e
# ╟─10f4a390-d6aa-4de9-a882-aa68e48064d9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
