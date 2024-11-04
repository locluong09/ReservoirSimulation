using DelimitedFiles
# using Plots
using LinearAlgebra

function calcTran_X(k1, k2, h1, h2, delta_x1, delta_x2, delta_y1, delta_y2, miu, B;beta_c = 1.127e-3)
    Ax1 = delta_y1*h1
    Ax2 = delta_y2*h2
    return Tran_X = 2*beta_c/(miu*B)/(delta_x1/(Ax1*k1) + delta_x2/(Ax2*k2))
end

function calcTran_Y(k1, k2, h1, h2, delta_x1, delta_x2, delta_y1, delta_y2, miu, B; beta_c = 1.127e-3)
    Ay1 = delta_x1*h1
    Ay2 = delta_x2*h2
    return 2*beta_c/(miu*B)/(delta_y1/(Ay1*k1) + (delta_y2/(Ay2*k2)))
end

function calcGamma(delta_x, delta_y, h; delta_t = 5, phi = 0.2, c = 3e-6, α_c = 5.615)
    return delta_x*delta_y*h*phi/(α_c*delta_t)
end

function calc_X(k1, k2, h1, h2, delta_x1, delta_x2, delta_y1, delta_y2; beta_c = 1.127e-3)
    #geometric factor
    Ax1 = delta_y1*h1
    Ax2 = delta_y2*h2
    return 2*beta_c/(delta_x1/(Ax1*k1) + delta_x2/(Ax2*k2))
end

function calc_Y(k1, k2, h1, h2, delta_x1, delta_x2, delta_y1, delta_y2; beta_c = 1.127e-3)
    Ay1 = delta_x1*h1
    Ay2 = delta_x2*h2
    return 2*beta_c/(delta_y1/(Ay1*k1) + delta_y2/(Ay2*k2))
end

kx, ky = 100, 100 # md
delta_x, delta_y = 500, 500 #ft
h = 50 #ft
ϕ = 0.2  #ft`
c_r = 3e-6 #psi-1
Pi = 4800 #psi
Swi = 0.3
Sgi = 0.1
Qo = 20 #STB/do
delta_t = 5 #day


oil_prop  = readdlm("/Users/luongkhanhloc/Desktop/CP3_Julia/oil_properties.txt")
water_prop  = readdlm("/Users/luongkhanhloc/Desktop/CP3_Julia/water_properties.txt")
gas_prop  = readdlm("/Users/luongkhanhloc/Desktop/CP3_Julia/gas_properties.txt")
oil_water_rel  = readdlm("/Users/luongkhanhloc/Desktop/CP3_Julia/oil_water_rel_table.txt")
gas_oil_rel  = readdlm("/Users/luongkhanhloc/Desktop/CP3_Julia/gas_oil_rel_table.txt")


function linearInterpolation(x::AbstractVector, y::AbstractVector, z::Float64)
    #=
    This function is not the same with matlat interp1, it just return the
    value of new value of fluid properties when get the value of z
    z is a new value of pressure that need to be interpolated to find the values
    of fluid properties
    =#

    n = length(x)
    if z in x
        return y[findall(x-> x==z, x)[1]]
    end
    # finding the interval that contains the value of z
    left = 1
    right = n
    while left < right
        if z < (x[left] + x[right]) / 2
            global right -= 1
        else
            global left += 1
        end
    end
    # update left and right because after the while loop right = left
    if z > x[right]
        global right += 1
    else
        global left -= 1
    end
    # using linear interpolation to update the value of z
    return (z - x[left])*(y[right] - y[left])/(x[right] - x[left]) + y[left]
end


function DlinearInterpolation(x::AbstractVector, y::AbstractVector, z::Float64)
    #=
    This function is not the same with matlat interp1, it just return the
    value of new value of fluid properties when get the value of z
    z is a new value of pressure that need to be interpolated to find the values
    of fluid properties
    =#

    n = length(x)
    if z in x
        left = findall(x-> x==z, x)[1]
        right = left + 1
        return (y[right] - y[left])/(x[right] - x[left])
    end
    # finding the interval that contains the value of z
    left = 1
    right = n
    while left < right
        if z < (x[left] + x[right]) / 2
            global right -= 1
        else
            global left += 1
        end
    end
    # update left and right because after the while loop right = left
    if z > x[right]
        global right += 1
    else
        global left -= 1
    end
    # using linear interpolation to update the value of z
    return (y[right] - y[left])/(x[right] - x[left])
end

function boundaries(X::AbstractArray)
    row, col = size(X)
    X_new  = zeros(row+2, col +2)
    for i in 1:row
        for j in 1:col
            X_new[i+1,j+1] = X[i,j]
        end
    end
    return X_new
end

Perm = ones(3,3)*100
Thick = ones(3,3)*50
DX = ones(3,3)*500
DY = ones(3,3)*500

Perm = boundaries(Perm)
Thickness = boundaries(Thick)
DX = boundaries(DX)
DY = boundaries(DY)
Row, Col = size(Perm)
cnt = 0
Order = zeros(Row, Col)
for i in 1:Row
    for j in 1:Col
        if Perm[i,j] > 0
            global cnt += 1
            Order[i,j] = cnt
        end
    end
end

Order = convert(Array{Int64}, Order)


Pre_init = ones(3,3)*4800
Sw_init = ones(3,3)*0.3
Sg_init = ones(3,3)*0.1

Pre_init = boundaries(Pre_init)
Sw_init = boundaries(Sw_init)
Sg_init = boundaries(Sg_init)
#initial pressure of oil phase
P_oil = boundaries(ones(3,3)*4800)
Sw = boundaries(ones(3,3)*0.3)
Sg = boundaries(ones(3,3)*0.1)

k_ro_Swc = 1

#Construct the transmissibility matrix
function construct_Residual_matrix(P_oil, Sg, Sw, Order, P_init, Sg_init, Sw_init)
    Row, Col = size(Order)
    m = Row - 2
    n = Col - 2
    T = zeros((Row-2)*(Col-2)*3, (Row-2)*(Col-2)*3)
    Q = zeros(m*n*3)
    for i = 1:Row
        for j = 1:Col
            if Order[i,j] > 0
                Gamma = calcGamma(500,500,50)
                S_w = Sw[i,j]
                S_g = Sg[i,j]

                P_cow = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,4], S_w)
                P_cgo = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,4], S_g)

                B_o_rhs = linearInterpolation(oil_prop[:,1], oil_prop[:,3], P_oil[i,j])
                B_g_rhs = linearInterpolation(gas_prop[:,1], gas_prop[:,3], P_oil[i,j] + P_cgo)
                B_w_rhs = linearInterpolation(water_prop[:,1], water_prop[:,3], P_oil[i,j] - P_cow)

                R_so_rhs = linearInterpolation(oil_prop[:,1], oil_prop[:,5], P_oil[i,j])

                if Order[i-1,j] > 0
                    P_o = (P_oil[i-1,j] + P_oil[i,j])/2
                    S_w_North = Sw[i-1,j]
                    S_g_North = Sg[i-1,j]
                    P_cow_North = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,4], S_w_North)
                    P_cgo_North = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,4], S_g_North)

                    #k_ro = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw)
                    if P_oil[i-1,j] > P_oil[i,j]
                        k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i-1,j])
                    else
                        k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i,j])
                    end

                    #k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw)
                    if P_oil[i-1,j] > P_oil[i,j]
                        k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i-1,j])
                    else
                        k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i,j])
                    end
                    #println(k_rw)
                    #k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg)
                    if P_oil[i-1,j] > P_oil[i,j]
                        k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i-1,j])
                    else
                        k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i,j])
                    end

                    #k_rog = inearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg)
                    if P_oil[i-1,j] > P_oil[i,j]
                        k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i-1,j])
                    else
                        k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i,j])
                    end
                    k_ro = k_ro_Swc*((k_row/k_ro_Swc + k_rw)*(k_rog/k_ro_Swc + k_rg) - (k_rw + k_rg))
                    B_o = linearInterpolation(oil_prop[:,1], oil_prop[:,3], P_o)
                    μ_o = linearInterpolation(oil_prop[:,1], oil_prop[:,4], P_o)

                    B_w = linearInterpolation(water_prop[:,1], water_prop[:,3], P_o - (P_cow+P_cow_North)/2)
                    μ_w = linearInterpolation(water_prop[:,1], water_prop[:,4], P_o - (P_cow+P_cow_North)/2)

                    B_g = linearInterpolation(gas_prop[:,1], gas_prop[:,3], P_o + (P_cgo+P_cgo_North)/2)
                    μ_g = linearInterpolation(gas_prop[:,1], gas_prop[:,4], P_o + (P_cgo+P_cgo_North)/2)
                    #first row is oil
                    North_o = calcTran_X(Perm[i-1,j], Perm[i,j], Thickness[i-1,j], Thickness[i,j], DX[i-1,j], DX[i,j], DY[i-1,j], DY[i,j], μ_o, B_o)
                    T[Order[i,j]*3-2, (Order[i-1,j]-1)*3 + 1] = North_o*k_ro
                    #second row is gas
                    North_g = calcTran_X(Perm[i-1,j], Perm[i,j], Thickness[i-1,j], Thickness[i,j], DX[i-1,j], DX[i,j], DY[i-1,j], DY[i,j], μ_g, B_g)
                    T[Order[i,j]*3-1, (Order[i-1,j]-1)*3 + 1] += North_g*k_rg
                    #third row is water
                    North_w = calcTran_X(Perm[i-1,j], Perm[i,j], Thickness[i-1,j], Thickness[i,j], DX[i-1,j], DX[i,j], DY[i-1,j], DY[i,j], μ_w, B_w)
                    T[Order[i,j]*3, (Order[i-1,j]-1)*3 + 1] = North_w*k_rw


                    #P_o Sg  and Sw component in matrix of coefficient in 1st row
                    T[Order[i,j]*3-2, Order[i,j]*3-2] -= North_o*k_ro


                    #P_o Sg  and Sw component in matrix of coefficient in 2nd row
                    R_so_lhs = linearInterpolation(oil_prop[:,1], oil_prop[:,5], P_o)
                    North_g_rso = North_o*k_ro*R_so_lhs
                    T[Order[i,j]*3-1, (Order[i-1,j]-1)*3 + 1] += North_g_rso
                    T[Order[i,j]*3-1, Order[i,j]*3-2] -= (North_g*k_rg + North_g_rso)


                    #P_o Sg  and Sw component in matrix of coefficient in 3rd row
                    T[Order[i,j]*3, Order[i,j]*3-2] -= North_w*k_rw
                    #gas is missing

                    Q[Order[i,j]*3-1] += (-P_cgo_North + P_cgo)*North_g*k_rg
                    Q[Order[i,j]*3-0] += (P_cow_North - P_cow)*North_w*k_rw


                end
                if Order[i+1,j] > 0
                    P_o = (P_oil[i+1,j] + P_oil[i,j])/2
                    S_w_South = Sw[i+1,j]
                    S_g_South = Sg[i+1,j]
                    P_cow_South = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,4], S_w_South)
                    P_cgo_South = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,4], S_g_South)


                    #k_ro = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw)
                    if P_oil[i+1,j] > P_oil[i,j]
                        k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i+1,j])
                    else
                        k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i,j])
                    end

                    #k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw)
                    if P_oil[i+1,j] > P_oil[i,j]
                        k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i+1,j])
                    else
                        k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i,j])
                    end

                    #k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg)
                    if P_oil[i+1,j] > P_oil[i,j]
                        k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i+1,j])
                    else
                        k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i,j])
                    end

                    #k_rog = inearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg)
                    if P_oil[i+1,j] > P_oil[i,j]
                        k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i+1,j])
                    else
                        k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i,j])
                    end
                    k_ro = k_ro_Swc*((k_row/k_ro_Swc + k_rw)*(k_rog/k_ro_Swc + k_rg) - (k_rw + k_rg))
                    B_o = linearInterpolation(oil_prop[:,1], oil_prop[:,3], P_o)
                    μ_o = linearInterpolation(oil_prop[:,1], oil_prop[:,4], P_o)

                    B_w = linearInterpolation(water_prop[:,1], water_prop[:,3], P_o - (P_cow+P_cow_South)/2)
                    μ_w = linearInterpolation(water_prop[:,1], water_prop[:,4], P_o - (P_cow+P_cow_South)/2)

                    B_g = linearInterpolation(gas_prop[:,1], gas_prop[:,3], P_o + (P_cgo+P_cgo_South)/2)
                    μ_g = linearInterpolation(gas_prop[:,1], gas_prop[:,4], P_o + (P_cgo+P_cgo_South)/2)
                    #first row is oil
                    South_o = calcTran_X(Perm[i+1,j], Perm[i,j], Thickness[i+1,j], Thickness[i,j], DX[i+1,j], DX[i,j], DY[i+1,j], DY[i,j], μ_o, B_o)
                    T[Order[i,j]*3-2, (Order[i+1,j]-1)*3 + 1] = South_o*k_ro
                    #second row is gas
                    South_g = calcTran_X(Perm[i+1,j], Perm[i,j], Thickness[i+1,j], Thickness[i,j], DX[i+1,j], DX[i,j], DY[i+1,j], DY[i,j], μ_g, B_g)
                    T[Order[i,j]*3-1, (Order[i+1,j]-1)*3 + 1] += South_g*k_rg
                    #third row is water
                    South_w = calcTran_X(Perm[i+1,j], Perm[i,j], Thickness[i+1,j], Thickness[i,j], DX[i+1,j], DX[i,j], DY[i+1,j], DY[i,j], μ_w, B_w)
                    T[Order[i,j]*3, (Order[i+1,j]-1)*3 + 1] = South_w*k_rw

                    #P_o Sg  and Sw component in matrix of coefficient in 1st row
                    T[Order[i,j]*3-2, Order[i,j]*3-2] -= South_o*k_ro


                    #P_o Sg  and Sw component in matrix of coefficient in 2nd row
                    R_so_lhs = linearInterpolation(oil_prop[:,1], oil_prop[:,5], P_o)
                    South_g_rso = South_o*k_ro*R_so_lhs
                    T[Order[i,j]*3-1, (Order[i+1,j]-1)*3 + 1] += South_g_rso
                    T[Order[i,j]*3-1, Order[i,j]*3-2] -= (South_g*k_rg + South_g_rso)


                    #P_o Sg  and Sw component in matrix of coefficient in 3rd row
                    T[Order[i,j]*3, Order[i,j]*3-2] -= South_w*k_rw
                    #gas is missing
                    Q[Order[i,j]*3-1] += (-P_cgo_South + P_cgo)*South_g*k_rg
                    Q[Order[i,j]*3-0] += (P_cow_South - P_cow)*South_w*k_rw


                end
                if Order[i,j-1] > 0
                    P_o = (P_oil[i,j-1] + P_oil[i,j])/2
                    S_w_West = Sw[i,j-1]
                    S_g_West = Sg[i,j-1]
                    P_cow_West = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,4], S_w_West)
                    P_cgo_West = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,4], S_g_West)


                    #k_ro = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw)
                    if P_oil[i,j-1] > P_oil[i,j]
                        k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i,j-1])
                    else
                        k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i,j])
                    end

                    #k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw)
                    if P_oil[i,j-1] > P_oil[i,j]
                        k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i,j-1])
                    else
                        k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i,j])
                    end

                    #k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg)
                    if P_oil[i,j-1] > P_oil[i,j]
                        k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i,j-1])
                    else
                        k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i,j])
                    end

                    #k_rog = inearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg)
                    if P_oil[i,j-1] > P_oil[i,j]
                        k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i,j-1])
                    else
                        k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i,j])
                    end
                    k_ro = k_ro_Swc*((k_row/k_ro_Swc + k_rw)*(k_rog/k_ro_Swc + k_rg) - (k_rw + k_rg))
                    B_o = linearInterpolation(oil_prop[:,1], oil_prop[:,3], P_o)
                    μ_o = linearInterpolation(oil_prop[:,1], oil_prop[:,4], P_o)

                    B_w = linearInterpolation(water_prop[:,1], water_prop[:,3], P_o - (P_cow+P_cow_West)/2)
                    μ_w = linearInterpolation(water_prop[:,1], water_prop[:,4], P_o - (P_cow+P_cow_West)/2)

                    B_g = linearInterpolation(gas_prop[:,1], gas_prop[:,3], P_o + (P_cgo+P_cgo_West)/2)
                    μ_g = linearInterpolation(gas_prop[:,1], gas_prop[:,4], P_o + (P_cgo+P_cgo_West)/2)
                    #first row is oil
                    West_o = calcTran_X(Perm[i,j-1], Perm[i,j], Thickness[i,j-1], Thickness[i,j], DX[i,j-1], DX[i,j], DY[i,j-1], DY[i,j], μ_o, B_o)
                    T[Order[i,j]*3-2, (Order[i,j-1]-1)*3 + 1] = West_o*k_ro
                    #second row is gas
                    West_g = calcTran_X(Perm[i,j-1], Perm[i,j], Thickness[i,j-1], Thickness[i,j], DX[i,j-1], DX[i,j], DY[i,j-1], DY[i,j], μ_g, B_g)
                    T[Order[i,j]*3-1, (Order[i,j-1]-1)*3 + 1] += West_g*k_rg
                    #third row is water
                    West_w = calcTran_X(Perm[i,j-1], Perm[i,j], Thickness[i,j-1], Thickness[i,j], DX[i,j-1], DX[i,j], DY[i,j-1], DY[i,j], μ_w, B_w)
                    T[Order[i,j]*3, (Order[i,j-1]-1)*3 + 1] = West_w*k_rw

                    #P_o Sg  and Sw component in matrix of coefficient in 1st row
                    T[Order[i,j]*3-2, Order[i,j]*3-2] -= West_o*k_ro

                    #P_o Sg  and Sw component in matrix of coefficient in 2nd row
                    R_so_lhs = linearInterpolation(oil_prop[:,1], oil_prop[:,5], P_o)
                    West_g_rso = West_o*k_ro*R_so_lhs
                    T[Order[i,j]*3-1, (Order[i,j-1]-1)*3 + 1] += West_g_rso
                    T[Order[i,j]*3-1, Order[i,j]*3-2] -= (West_g*k_rg + West_g_rso)


                    #P_o Sg  and Sw component in matrix of coefficient in 3rd row
                    T[Order[i,j]*3, Order[i,j]*3-2] -= West_w*k_rw

                    Q[Order[i,j]*3-1] += (-P_cgo_West + P_cgo)*West_g*k_rg
                    Q[Order[i,j]*3-0] += (P_cow_West - P_cow)*West_w*k_rw

                end
                if Order[i,j+1] > 0
                    P_o = (P_oil[i,j+1] + P_oil[i,j])/2
                    S_w_East = Sw[i,j+1]
                    S_g_East = Sg[i,j+1]
                    P_cow_East = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,4], S_w_East)
                    P_cgo_East = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,4], S_g_East)


                    #k_ro = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw)
                    if P_oil[i,j+1] > P_oil[i,j]
                        k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i,j+1])
                    else
                        k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i,j])
                    end

                    #k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw)
                    if P_oil[i,j+1] > P_oil[i,j]
                        k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i,j+1])
                    else
                        k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i,j])
                    end

                    #k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg)
                    if P_oil[i,j+1] > P_oil[i,j]
                        k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i,j+1])
                    else
                        k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i,j])
                    end

                    #k_rog = inearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg)
                    if P_oil[i,j+1] > P_oil[i,j]
                        k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i,j+1])
                    else
                        k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i,j])
                    end
                    k_ro = k_ro_Swc*((k_row/k_ro_Swc + k_rw)*(k_rog/k_ro_Swc + k_rg) - (k_rw + k_rg))
                    B_o = linearInterpolation(oil_prop[:,1], oil_prop[:,3], P_o)
                    μ_o = linearInterpolation(oil_prop[:,1], oil_prop[:,4], P_o)

                    B_w = linearInterpolation(water_prop[:,1], water_prop[:,3], P_o - (P_cow+P_cow_East)/2)
                    μ_w = linearInterpolation(water_prop[:,1], water_prop[:,4], P_o - (P_cow+P_cow_East)/2)

                    B_g = linearInterpolation(gas_prop[:,1], gas_prop[:,3], P_o + (P_cgo+P_cgo_East)/2)
                    μ_g = linearInterpolation(gas_prop[:,1], gas_prop[:,4], P_o + (P_cgo+P_cgo_East)/2)
                    #first row is oil
                    East_o = calcTran_X(Perm[i,j+1], Perm[i,j], Thickness[i,j+1], Thickness[i,j], DX[i,j+1], DX[i,j], DY[i,j+1], DY[i,j], μ_o, B_o)
                    T[Order[i,j]*3-2, (Order[i,j+1]-1)*3 + 1] = East_o*k_ro
                    #second row is gas
                    East_g = calcTran_X(Perm[i,j+1], Perm[i,j], Thickness[i,j+1], Thickness[i,j], DX[i,j+1], DX[i,j], DY[i,j+1], DY[i,j], μ_g, B_g)
                    T[Order[i,j]*3-1, (Order[i,j+1]-1)*3 + 1] += East_g*k_rg
                    #third row is water
                    East_w = calcTran_X(Perm[i,j+1], Perm[i,j], Thickness[i,j+1], Thickness[i,j], DX[i,j+1], DX[i,j], DY[i,j+1], DY[i,j], μ_w, B_w)
                    T[Order[i,j]*3, (Order[i,j+1]-1)*3 + 1] = East_w*k_rw

                    #P_o Sg  and Sw component in matrix of coefficient in 1st row
                    T[Order[i,j]*3-2, Order[i,j]*3-2] -= East_o*k_ro

                    #P_o Sg  and Sw component in matrix of coefficient in 2nd row
                    R_so_lhs = linearInterpolation(oil_prop[:,1], oil_prop[:,5], P_o)
                    East_g_rso = East_o*k_ro*R_so_lhs
                    T[Order[i,j]*3-1, (Order[i,j+1]-1)*3 + 1] += East_g_rso
                    T[Order[i,j]*3-1, Order[i,j]*3-2] -= (East_g*k_rg + East_g_rso)


                    #P_o Sg  and Sw component in matrix of coefficient in 3rd row
                    T[Order[i,j]*3, Order[i,j]*3-2] -= East_w*k_rw

                    Q[Order[i,j]*3-1] += (-P_cgo_East + P_cgo)*East_g*k_rg
                    Q[Order[i,j]*3-0] += (P_cow_East - P_cow)*East_w*k_rw

                end


                T[Order[i,j]*3-2, Order[i,j]*3-1] = Gamma/B_o_rhs
                T[Order[i,j]*3-2, Order[i,j]*3-0] = Gamma/B_o_rhs

                T[Order[i,j]*3-1, Order[i,j]*3-1] = Gamma*(R_so_rhs/B_o_rhs - 1/B_g_rhs)
                T[Order[i,j]*3-1, Order[i,j]*3-0] = Gamma*R_so_rhs/B_o_rhs

                T[Order[i,j]*3, Order[i,j]*3-0] = -Gamma/B_w_rhs
                #RHS of Q
                B_o_init = linearInterpolation(oil_prop[:,1], oil_prop[:,3], P_init[i,j])
                B_w_init = linearInterpolation(water_prop[:,1], water_prop[:,3], P_init[i,j])
                B_g_init = linearInterpolation(gas_prop[:,1], gas_prop[:,3], P_init[i,j])
                # println(B_g_init)
                R_so_init = linearInterpolation(oil_prop[:,1], oil_prop[:,5], P_init[i,j])

                # Q[Order[i,j]*3-2] = Gamma*(1/B_o_rhs - B_o_init) + Gamma*(Sw_init[i,j] + Sg_init[i,j])/B_o_init
                # Q[Order[i,j]*3-1] += R_so_rhs*(Gamma*(1/B_o_rhs - 1/B_o_init) + Gamma*(Sw_init[i,j] + Sg_init[i,j])/B_o_init) - Gamma*(Sg_init[i,j]/B_g_init)
                # Q[Order[i,j]*3-0] += -Gamma*(Sw_init[i,j]/B_w_init)
                Q[Order[i,j]*3-2] += Gamma*(1/B_o_rhs - 1/B_o_init + (Sw_init[i,j] + Sg_init[i,j])/B_o_init)
                Q[Order[i,j]*3-1] += Gamma*((R_so_rhs/B_o_rhs - R_so_init/B_o_init) + R_so_init*(Sw_init[i,j] + Sg_init[i,j])/B_o_init) - Gamma*(Sg_init[i,j]/B_g_init)
                Q[Order[i,j]*3-0] += -Gamma*(Sw_init[i,j]/B_w_init)

            end
        end
    end
    for i = 1:Row
        for j = 1:Col
            if Order[i,j] == 5
                re = 0.198*DX[i,j]
                rw = 0.25
                P_cow = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,4], Sw[i,j])
                P_cgo = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,4], Sg[i,j])

                B_o_well = linearInterpolation(oil_prop[:,1], oil_prop[:,3], P_oil[i,j])
                μ_o_well = linearInterpolation(oil_prop[:,1], oil_prop[:,4], P_oil[i,j])

                B_w_well = linearInterpolation(water_prop[:,1], water_prop[:,3], P_oil[i,j] - P_cow)
                μ_w_well = linearInterpolation(water_prop[:,1], water_prop[:,4], P_oil[i,j] - P_cow)

                B_g_well = linearInterpolation(gas_prop[:,1], gas_prop[:,3], P_oil[i,j] + P_cgo)
                μ_g_well = linearInterpolation(gas_prop[:,1], gas_prop[:,4], P_oil[i,j] + P_cgo)

                k_row = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,3], Sw[i,j])
                k_rw = linearInterpolation(oil_water_rel[:,1], oil_water_rel[:,2], Sw[i,j])
                k_rog = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,3], Sg[i,j])
                k_rg = linearInterpolation(gas_oil_rel[:,1], gas_oil_rel[:,2], Sg[i,j])
                R_so = linearInterpolation(oil_prop[:,1], oil_prop[:,5], P_oil[i,j])
                k_ro_well = k_ro = k_ro_Swc*(k_row/k_ro_Swc + k_rw)*(k_rog/k_ro_Swc + k_rg) - (k_rw + k_rg)
                P_sf = P_oil[i,j] - 20/(1.127*1e-3*2*pi*Perm[i,j]*50/log(re/rw))/(k_ro_well/B_o_well/μ_o_well)
                Q[Order[i,j]*3 - 2] += 20
                Q[Order[i,j]*3 - 1] += (1.127*1e-3*2*pi*Perm[i,j]*50/log(re/rw))*(k_rg/(μ_g_well*B_g_well) + R_so*k_ro_well/(μ_o_well*B_o_well))*(P_oil[i,j] - P_sf)+
                                        R_so*20
                Q[Order[i,j]*3 - 0] += (1.127*1e-3*2*pi*Perm[i,j]*50/log(re/rw))*(k_rw/(μ_w_well*B_w_well))*(P_oil[i,j] - P_sf)
            end
        end
    end

    return T,Q
end



function Jacobian(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init)
    Row, Col = size(Order)
    m = Row - 2
    n = Col - 2
    x = zeros(m*n*3)
    for i in 1:Row
        for j in 1:Col
            if Order[i,j] > 0
                x[Order[i,j]*3-2] = P_oil[i,j]
                x[Order[i,j]*3-1] = Sg[i,j]
                x[Order[i,j]*3] = Sw[i,j]
            end
        end
    end
    δp = 1e-3
    δSg = 1e-7
    δSw = 1e-7

    J = zeros(27,27)
    T,Q = construct_Residual_matrix(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init)
    Residual = T*x - Q
    for i in 1:Row
        for j in 1:Col
            if Order[i,j] > 0
                P_oil_new = copy(P_oil)
                P_oil_new[i,j] += δp
                #println(P_oil_new)

                x_new1 = copy(x)
                x_new1[Order[i,j]*3-2] += δp

                To,Qo = construct_Residual_matrix(P_oil_new, Sg, Sw, Order, Pre_init, Sg_init, Sw_init)
                R_new1 = To*x_new1 - Qo
                #println(R_new1 - Residual)
                # for k = 1:12
                #     J[k,Order[i,j]*3-2] = (R_new[k] - Residual[k])/δp
                # end
                J[:,Order[i,j]*3-2] = (R_new1 - Residual)/δp
                #println(J[:,Order[i,j]*3-2])

                Sg_new = copy(Sg)
                Sg_new[i,j] += δSg

                x_new2 = copy(x)
                x_new2[Order[i,j]*3-1] += δSg

                Tg,Qg = construct_Residual_matrix(P_oil, Sg_new, Sw, Order, Pre_init, Sg_init, Sw_init)
                R_new2 = Tg*x_new2 - Qg
                # for k = 1:12
                #     J[k,Order[i,j]*3-1] = (R_new[k] - Residual[k])/δSg
                # end
                J[:,Order[i,j]*3-1] = (R_new2 - Residual)/δSg


                Sw_new = copy(Sw)
                Sw_new[i,j] += δSw

                x_new3 = copy(x)
                x_new3[Order[i,j]*3] += δSw

                Tw,Qw = construct_Residual_matrix(P_oil, Sg, Sw_new, Order, Pre_init, Sg_init, Sw_init)
                R_new3 = Tw*x_new3 - Qw
                # for k = 1:12
                #     J[k,Order[i,j]*3] = (R_new[k] - Residual[k])/δSw
                # end
                J[:,Order[i,j]*3] = (R_new3 - Residual)/δSw
            end


        end
    end

    return J
end


x = zeros(27,1)

for i in 1:Row
    for j in 1:Col
        if Order[i,j] > 0
            x[Order[i,j]*3-2] = P_oil[i,j]
            x[Order[i,j]*3-1] = Sg[i,j]
            x[Order[i,j]*3] = Sw[i,j]
        end
    end
end
iter = 0

for i in 1:3
# while true
    #calculate everything for the loop
    #Initiaze the pressure again
    global iter += 1
    # println(x)
    T, Q  = construct_Residual_matrix(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init)
    Residual = T*x - Q
    # println(Residual)
    J = Jacobian(P_oil, Sg, Sw, Order, Pre_init, Sg_init, Sw_init)
    println("Jacobian matrix in iteration ", i)
    println(J)
    
    Delta = J\Residual


    x_new = x
    
    x_new = x_new - Delta
    
    P_update = zeros(Row,Col)
    Sg_update = zeros(Row,Col)
    Sw_update = zeros(Row,Col)
    #print(x_new)
    for i in 1:Row
        for j in 1:Col
            if Order[i,j] > 0
                P_update[i,j] = x_new[Order[i,j]*3-2]
                Sg_update[i,j] = x_new[Order[i,j]*3-1]
                Sw_update[i,j] = x_new[Order[i,j]*3]
            end
        end
    end
    # println(maximum(broadcast(abs, x_new-x)))
    # if maximum(broadcast(abs,Delta)) < 0.01
    #     break
    # end
    # if maximum(broadcast(abs, x_new-x)) < 0.01
    #     break
    # end

    global P_oil = P_update
    global Sg = Sg_update
    global Sw = Sw_update
    #Update the pressure
    #Pre_init, Pre_predict = Pre_predict, Pre_update
    #Update Z
    #Z_predict_local = gasCompressibility(Pre_update, 150., Order)
    #global Z_predict = Z_predict_local

    #println(maximum(Pre_predict))
    #println(Pre_predict)
    global x = x_new

end

println(iter)
