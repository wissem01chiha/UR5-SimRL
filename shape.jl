"""
 shape_function : return lagrange shape functions

 Args   :  
    - X1 (float) local elment referance 
    - X2 (float) local elment refrenace 
    - K (integer) shape function index number
    - N (integer) x coordinate vector length

 Author :
    wissem chiha   
 Date   :
    February 2024 
"""
function shape_function(x1::double, x2::double, k::Int, n::int)::Vector{float}
    he = x2 - x1
    x = range(0, he, length=n)
    l = similar(x)
    if k == 1
        l .= 1 .- 3 .* (x.^2) ./ (he^2) + 2 .* (x.^3) ./ he^3
    elseif k == 2
        l .= x .- 2 .* (x.^2) ./ he .+ (x.^3) ./ (he^2)
    elseif k == 3
        l .= 3 .* (x.^2) ./ (he^2) .- 2 .* (x.^3) ./ he^3
    elseif k == 4
        l .= -(x.^2) ./ he .+ (x.^3) ./ (he^2)
    end 
    
    return l
end

"""
 d2_shape_function: return second derivitative lagrange 
                    shape functions
 Args:  
    - X1 (float) local elment referance 
    - X2 (float) local elment refrenace 
    - K (integer) shape function index number
    - N (integer) x coordinate vector length
 Date:
 

"""
function d2_shape_function(x1::double, x2::double, k::int, n::int)::Vector{float}
    he = x2 - x1
    x = range(0, he, length=n)
    l = l = similar(x)
    if k == 1
        l = -6 / (he^2) .+ 12 .* x ./ he^3 
    elseif k == 2
        l = -4 / he .+ 6 .* (x ./ he.^2)
    elseif k == 3
        l = 6 / (he^2) .- 12 .* (x ./ he^3)
    elseif k == 4
        l = -2 / he .+ 6 .* (x ./ he^2)
    end
    return l
end



