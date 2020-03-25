using Plots
using DelimitedFiles
using PyPlot
using Plotly
#include("boundaries.jl")

perm  = readdlm("/Users/luongkhanhloc/Desktop/CP2-Julia/k.txt")
top  = readdlm("/Users/luongkhanhloc/Desktop/CP2-Julia/t.txt")
thickness  = readdlm("/Users/luongkhanhloc/Desktop/CP2-Julia/h.txt")
dx  = readdlm("/Users/luongkhanhloc/Desktop/CP2-Julia/dx.txt")
dy  = readdlm("/Users/luongkhanhloc/Desktop/CP2-Julia/dy.txt")
x_coordinate  = readdlm("/Users/luongkhanhloc/Desktop/CP2-Julia/x_coordinate.txt")
y_coordinate  = readdlm("/Users/luongkhanhloc/Desktop/CP2-Julia/y_coordinate.txt")

function meshgrid(x,y)
    m = length(x)
    n = length(y)
    vx = reshape(x, 1 , m)
    vy = reshape(y, n ,1)
    return repeat(vx,n,1), repeat(vy,1,m)
end
X,Y = meshgrid(x_coordinate, y_coordinate)
surf(X,Y,perm)

#Change all negative values to be zero
row, col = size(perm)
for i = 1:row
    for j = 1:col
        if perm[i,j] < 0
            perm[i,j] = 0
        end
        if thickness[i,j] < 0
            thickness[i,j] = 0
        end
        if top[i,j] > 3100
                top[i,j] = 3100
        end
    end
end
order = zeros(Int64, (row,col))
cnt = 0
for i in 1:row
    for j in 1:col
        if perm[i,j] > 0
            global cnt += 1
            order[i,j] = cnt
        end
    end
end


function create_matrix_delta(X::AbstractArray, Y::AbstractArray)
    m = length(X)
    n = length(Y)
    DX = zeros(m-1)
    DY = zeros(n-1)
    for i = 1:length(DX)
        DX[i] = X[i+1] - X[i]
    end
    for i = 1:length(DY)
        DY[i] = Y[i+1] - Y[i]
    end
    vx = reshape(DX,1,m-1)
    vy = reshape(DY,n-1,1)
    repeat(vx,n-1,1), repeat(vy,1,m-1)
end
DX,DY = create_matrix_delta(dx,dy)



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

Perm = boundaries(perm)
Top = boundaries(top)
Thickness = boundaries(thickness)
Order = boundaries(order)
Order = convert(Array{Int64}, Order)
DX = boundaries(DX)
DY = boundaries(DY)
Row,Col = size(Perm)

for i in 1:Row
    for j = 1:Col
        if Order[i,j] == 6
            Thickness[i,j] = (Thickness[i-1,j] + Thickness[i+1,j] + Thickness[i,j-1] + Thickness[i,j+1])/4
        end
    end
end

#Compressibility factor from Dranchuk and Abou-Kassem, 1975
function gasCompressibility(P,index; T = 150, Ppc = 680., Tpc = 385.)
    Pr = P/Ppc
    Tr = (T+460)/Tpc
	A=[.3265; -1.07; -0.5339; 0.01569; -0.05165; 0.5475; -0.7361; 0.1844; 0.1056;
	 0.6134; 0.721]
	c1 = A[1] + A[2]/Tr + A[3]*Tr^-3 + A[4]*Tr^-4 + A[5]*Tr^-5
	c2 = A[6] + A[7]/Tr + A[8]/Tr^2
	c3 = A[7]/Tr + A[8]/Tr^2
	u,v = size(index,1), size(index,2)
	ro = zeros(u,v)
	ron = zeros(u,v)
	z = zeros(u,v)
	for i = 1:u
		for j = 1:v
			if index[i,j] != 0
				ro[i,j] = 0.27*Pr[i,j]/Tr
				diff = 1
				while diff > 0.01
					f = 1-0.27*Pr[i,j]/(ro[i,j]*Tr) + c1*ro[i,j] +
                    c2*ro[i,j]^2 - A[9]*c3*ro[i,j]^5 +
                    A[10]*(1 + A[11]*ro[i,j]^2)*ro[i,j]^2*exp(-1*A[11]*ro[i,j]^2)/Tr^3;

                    fd = 0.27*Pr[i,j]/(ro[i,j]^2*Tr) + c1 + 2*c2*ro[i,j] -
                    5*A[9]*c3*ro[i,j]^4 +
                    2*A[10]*ro[i,j]*(1+A[11]*ro[i,j]^2 -
                    A[11]^2*ro[i,j]^4)*exp(-1*A[11]*ro[i,j]^2)/Tr^3;

                    ron[i,j] = ro[i,j] - f/fd
                    diff = abs(ron[i,j] - ro[i,j])
                    ro[i,j] = ron[i,j]
                end
                z[i,j] = 0.27*Pr[i,j]/(ro[i,j]*Tr)
            end
        end
    end
    if length(z) == 1
        return z[1]
    else
        return z
    end
end

function DgasCompressibility(P, index; T = 150, Ppc = 680., Tpc = 385.)
    Pr = P/Ppc
    Tr = (T+460)/Tpc
	A=[.3265; -1.07; -0.5339; 0.01569; -0.05165; 0.5475; -0.7361; 0.1844; 0.1056;
	 0.6134; 0.721]
	c1 = A[1] + A[2]/Tr + A[3]*Tr^-3 + A[4]*Tr^-4 + A[5]*Tr^-5
	c2 = A[6] + A[7]/Tr + A[8]/Tr^2
	c3 = A[7]/Tr + A[8]/Tr^2
	u,v = size(index,1), size(index,2)
	ro = zeros(u,v)
	ron = zeros(u,v)
	z = zeros(u,v)
	for i = 1:u
		for j = 1:v
			if index[i,j] != 0
				ro[i,j] = 0.27*Pr[i,j]/Tr
				diff = 1
				while diff > 0.01
					f = 1-0.27*Pr[i,j]/(ro[i,j]*Tr) + c1*ro[i,j] +
                    c2*ro[i,j]^2 - A[9]*c3*ro[i,j]^5 +
                    A[10]*(1 + A[11]*ro[i,j]^2)*ro[i,j]^2*exp(-1*A[11]*ro[i,j]^2)/Tr^3;

                    fd = 0.27*Pr[i,j]/(ro[i,j]^2*Tr) + c1 + 2*c2*ro[i,j] -
                    5*A[9]*c3*ro[i,j]^4 +
                    2*A[10]*ro[i,j]*(1+A[11]*ro[i,j]^2 -
                    A[11]^2*ro[i,j]^4)*exp(-1*A[11]*ro[i,j]^2)/Tr^3;

                    ron[i,j] = ro[i,j] - f/fd
                    diff = abs(ron[i,j] - ro[i,j])
                    ro[i,j] = ron[i,j]
                end
                z[i,j] = 0.27*Pr[i,j]/(ro[i,j]*Tr)
            end
        end
    end
    if length(z) == 1
        return z[1]/P
    else
        return z/P
    end
end

function gasZ(p)
    z = 4.949892*1e-8*p^2 - 2.30415*1e-4p + 1.00563
    return z
end

function DgasZ(p)
    dz = 2*4.949892*1e-8*p - 2.30415*1e-4
    return dz
end


function gasDensity(pressure; SG= 0.6, R = 10.731, T = 150)
    T_new = 150+460
    MW = 0.6*28.96
    return ρ = pressure*MW/(R*T_new)
end

function DgasDensity(pressure; SG= 0.6, R = 10.731, T = 150)
    T_new = 150+460
    MW = 0.6*28.96
    return dρ = MW/(R*T_new)
end


#Gas viscosity function from Lee, Gonzalez and Eakin
function gasVis(ρ; T = 150., SG = 0.6, R = 10.731)
    MW = SG*28.96
    T_new = T + 460
    X = 3.448 +986.4/T_new + 0.01009*MW
    Y = 2.447 - 0.224X
    K = (9.379 + 0.01607MW)*T_new^1.5/(209.2+19.26MW +T_new)
    return μ = 1e-4*K*exp(X*(ρ/62.4)^Y)
end

function DgasVis(ρ; T = 150., SG = 0.6, R = 10.731)
    MW = SG*28.96
    T_new = T + 460
    X = 3.448 +986.4/T_new + 0.01009*MW
    Y = 2.447 - 0.224X
    K = (9.379 + 0.01607MW)*T_new^1.5/(209.2+19.26MW +T_new)
    return dμ = 1e-4*K*X*Y*(ρ/62.4)^(Y-1)*exp(X*(ρ/62.4)^Y)/62.4
end
#Gas formation volume factor
function gasFVF(pressure, z; T = 150, p_sc = 14.7, T_sc = 520)
    T_new = 150+460
    return B_g = p_sc/T_sc/5.615*z*T_new/pressure
end

function DgasFVF(pressure, z; T = 150, p_sc = 14.7, T_sc = 520)
    T_new = 150+460
    return DB_g = -p_sc/T_sc/5.615*T_new*(DgasZ(pressure)*pressure - gasZ(pressure))/pressure^2
end

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

function calcGamma(delta_x, delta_y, h; T = 150., delta_t = 30, phi = 1,p_sc = 14.7, T_sc = 520)
    return delta_x*delta_y*h*phi*T_sc/(p_sc*(T+460)*delta_t)
end

function calc_X(k1, k2, h1, h2, delta_x1, delta_x2, delta_y1, delta_y2; beta_c = 1.127e-3)
    Ax1 = delta_y1*h1
    Ax2 = delta_y2*h2
    return 2*beta_c/(delta_x1/(Ax1*k1) + delta_x2/(Ax2*k2))
end

function calc_Y(k1, k2, h1, h2, delta_x1, delta_x2, delta_y1, delta_y2; beta_c = 1.127e-3)
    Ay1 = delta_x1*h1
    Ay2 = delta_x2*h2
    return 2*beta_c/(delta_y1/(Ay1*k1) + delta_y2/(Ay2*k2))
end




Pre_init = zeros(Row,Col)
Pre_predict = zeros(Row, Col)
for i in 1:Row
    for j in 1:Col
        if Order[i,j] > 0
            Pre_init[i,j] = 4000
            Pre_predict[i,j] = 3000
        end
    end
end

Z_init = gasCompressibility(Pre_init, Order)
Z_predict = gasCompressibility(Pre_predict, Order)
Gamma = zeros(Row, Col)
for i in 1:Row
    for j in 1:Col
        if Order[i,j] > 0
            Gamma[i,j] = calcGamma(DX[i,j], DY[i,j], Thickness[i,j])
        end
    end
end


#Construct the transmissibility matrix
function construct_Tran_matrix(Pre_init, Pre_predict, Gamma, Oder_matrix)
    Row, Col = size(Order)
    N = zeros(Row, Col)
    S = zeros(Row, Col)
    W = zeros(Row, Col)
    E = zeros(Row, Col)
    C = zeros(Row, Col)
    Q = zeros(Row, Col)
    for i = 1:Row
        for j = 1:Col
            if Order[i,j] > 0
                if Order[i-1,j] > 0
                    Pressure = (Pre_predict[i-1,j] + Pre_predict[i,j])/2
                    ρ = gasDensity(Pressure)
                    μ = gasVis(ρ)
                    z = gasZ(Pressure)
                    Bg = gasFVF(Pressure, z)
                    N[i,j] = calcTran_X(Perm[i-1,j], Perm[i,j], Thickness[i-1,j], Thickness[i,j], DX[i-1,j], DX[i,j], DY[i-1,j], DY[i,j], μ, Bg)
                end
                if Order[i+1,j] > 0
                    Pressure = (Pre_predict[i+1,j] + Pre_predict[i,j])/2
                    ρ = gasDensity(Pressure)
                    μ = gasVis(ρ)
                    z = gasZ(Pressure)
                    Bg = gasFVF(Pressure, z)
                    S[i,j] = calcTran_X(Perm[i+1,j], Perm[i,j], Thickness[i+1,j], Thickness[i,j], DX[i+1,j], DX[i,j], DY[i+1,j], DY[i,j], μ, Bg)
                end
                if Order[i,j-1] > 0
                    Pressure = (Pre_predict[i,j-1] + Pre_predict[i,j])/2
                    ρ = gasDensity(Pressure)
                    μ = gasVis(ρ)
                    z = gasZ(Pressure)
                    Bg = gasFVF(Pressure, z)
                    W[i,j] = calcTran_Y(Perm[i,j-1], Perm[i,j], Thickness[i,j-1], Thickness[i,j], DX[i,j-1], DX[i,j], DY[i,j-1], DY[i,j], μ, Bg)
                end
                if Order[i,j+1] > 0
                    Pressure = (Pre_predict[i,j+1] + Pre_predict[i,j])/2
                    ρ = gasDensity(Pressure)
                    μ = gasVis(ρ)
                    z = gasZ(Pressure)
                    Bg = gasFVF(Pressure, z)
                    E[i,j] = calcTran_Y(Perm[i,j+1], Perm[i,j], Thickness[i,j+1], Thickness[i,j], DX[i,j+1], DX[i,j], DY[i,j+1], DY[i,j], μ, Bg)
                end
                C[i,j] = -(N[i,j] + S[i,j] + W[i,j] + E[i,j] + Gamma[i,j]/gasZ(Pre_predict[i,j]))
                Q[i,j] = Gamma[i,j]*(-Pre_init[i,j]/gasZ(Pre_init[i,j]))
            end
        end
    end
    Q[9,11] += 300e6
    Q[7,15] += 300e6
    return (N,S,W,E,C,Q)
end

#N,S,W,E,C,Q = construct_Tran_matrix(Pre_init, Pre_predict, Gamma, Z_init, Z_predict, Order)
#Calculate the Residual matrix
function construct_Coefficient_matrix(N,S,W,E,C,Q,Order,cnt)
    LHS = zeros(cnt, cnt)
    RHS = zeros(cnt, 1)
    for i = 1:Row
        for j = 1:Col
            if Order[i,j] > 0
                LHS[Order[i,j], Order[i,j]] = C[i,j]
                if Order[i-1,j] > 0
                    LHS[Order[i,j], Order[i-1,j]] = N[i,j]
                end
                if Order[i+1,j] > 0
                    LHS[Order[i,j], Order[i+1,j]] = S[i,j]
                end
                if Order[i,j-1] > 0
                    LHS[Order[i,j], Order[i,j-1]] = W[i,j]
                end
                if Order[i,j+1] > 0
                    LHS[Order[i,j], Order[i,j+1]] = E[i,j]
                end
            end
        end
    end
    for i in 1:Row
        for j in 1:Col
            if Order[i,j] > 0
                RHS[Order[i,j]] = Q[i,j]
            end
        end
    end
    return (LHS, RHS)
end

#LHS, RHS  = construct_Coefficient_matrix(N,S,W,E,C,Q,Order,cnt)
# Construct residual matrix
function Residual(LHS, RHS, Pre_predict, Order,cnt)
    Pre_var = zeros(cnt,1)
    Residual = zeros(cnt,1)
    for i in 1:Row
        for j in 1:Col
            if Order[i,j] > 0
                Pre_var[Order[i,j], 1] = Pre_predict[i,j]
            end
        end
    end
    # Residual matrix
    Residual = LHS*Pre_var- RHS
    return Residual
end

#Residual_matrix = Residual(LHS, RHS, Pre_predict, Order, cnt)
#Construct the Jacobian matrix
function Jacobian_matrix(Pressure, Order, Gamma, cnt)
    J = zeros(cnt,cnt)
    for i in 1:Row
        for j in 1:Col
            if Order[i,j] > 0
                if Order[i-1,j] > 0
                    N = calc_X(Perm[i-1,j], Perm[i,j], Thickness[i-1,j], Thickness[i,j], DX[i-1,j], DX[i,j], DY[i-1,j], DY[i,j])
                    Pressure = (Pre_predict[i-1,j] + Pre_predict[i,j])/2
                    ρ = gasDensity(Pressure)
                    μ = gasVis(ρ)
                    z = gasZ(Pressure)
                    Bg = gasFVF(Pressure, z)
                    Dμ = DgasVis(ρ)
                    Dρ = DgasDensity(Pressure)
                    DBg = DgasFVF(Pressure, z)
                    JN = 1/2*Pre_predict[i,j]*(Dμ*Dρ*Bg + μ*DBg)/(μ*Bg)^2 + (μ*Bg - Pre_predict[i-1,j]/2*(Dμ*Dρ*Bg + μ*DBg))/(μ*Bg)^2
                    J_N = JN * N
                    C_N = (μ*Bg - Pre_predict[i,j]/2*(Dμ*Dρ*Bg + μ*DBg))/(μ*Bg)^2*N
                    C_N1 = -1/2*Pre_predict[i-1,j]*(Dμ*Dρ*Bg + μ*DBg)/(μ*Bg)^2*N
                    J[Order[i,j], Order[i-1,j]] = J_N
                    J[Order[i,j], Order[i,j]] -= C_N
                    J[Order[i,j], Order[i,j]] += C_N1

                end
                if Order[i+1,j] > 0
                    S = calc_X(Perm[i+1,j], Perm[i,j], Thickness[i+1,j], Thickness[i,j], DX[i+1,j], DX[i,j], DY[i+1,j], DY[i,j])
                    Pressure = (Pre_predict[i+1,j] + Pre_predict[i,j])/2
                    ρ = gasDensity(Pressure)
                    μ = gasVis(ρ)
                    z = gasZ(Pressure)
                    Bg = gasFVF(Pressure, z)
                    Dμ = DgasVis(ρ)
                    Dρ = DgasDensity(Pressure)
                    DBg = DgasFVF(Pressure, z)
                    JS = 1/2*Pre_predict[i,j]*(Dμ*Dρ*Bg + μ*DBg)/(μ*Bg)^2 + (μ*Bg - Pre_predict[i+1,j]/2*(Dμ*Dρ*Bg + μ*DBg))/(μ*Bg)^2
                    J_S = JS * S
                    C_S = (μ*Bg - Pre_predict[i,j]/2*(Dμ*Dρ*Bg + μ*DBg))/(μ*Bg)^2*S
                    C_S1 = -1/2*Pre_predict[i+1,j]*(Dμ*Dρ*Bg + μ*DBg)/(μ*Bg)^2*S
                    J[Order[i,j], Order[i+1,j]] = J_S
                    J[Order[i,j], Order[i,j]] -= C_S
                    J[Order[i,j], Order[i,j]] += C_S1
                end
                if Order[i,j-1] > 0
                    W = calc_Y(Perm[i,j-1], Perm[i,j], Thickness[i,j-1], Thickness[i,j], DX[i,j-1], DX[i,j], DY[i,j-1], DY[i,j])
                    Pressure = (Pre_predict[i,j-1] + Pre_predict[i,j])/2
                    ρ = gasDensity(Pressure)
                    μ = gasVis(ρ)
                    z = gasZ(Pressure)
                    Bg = gasFVF(Pressure, z)
                    Dμ = DgasVis(ρ)
                    Dρ = DgasDensity(Pressure)
                    DBg = DgasFVF(Pressure, z)
                    JW = 1/2*Pre_predict[i,j]*(Dμ*Dρ*Bg + μ*DBg)/(μ*Bg)^2 + (μ*Bg - Pre_predict[i,j-1]/2*(Dμ*Dρ*Bg + μ*DBg))/(μ*Bg)^2
                    J_W = JW * W
                    C_W = (μ*Bg - Pre_predict[i,j]/2*(Dμ*Dρ*Bg + μ*DBg))/(μ*Bg)^2*W
                    C_W1 = -1/2*Pre_predict[i,j-1]*(Dμ*Dρ*Bg + μ*DBg)/(μ*Bg)^2*W
                    J[Order[i,j], Order[i,j-1]] = J_W
                    J[Order[i,j], Order[i,j]] -= C_W
                    J[Order[i,j], Order[i,j]] += C_W1
                end
                if Order[i,j+1] > 0
                    E = calc_Y(Perm[i,j+1], Perm[i,j], Thickness[i,j+1], Thickness[i,j], DX[i,j+1], DX[i,j], DY[i,j+1], DY[i,j])
                    Pressure = (Pre_predict[i,j+1] + Pre_predict[i,j])/2
                    ρ = gasDensity(Pressure)
                    μ = gasVis(ρ)
                    z = gasZ(Pressure)
                    Bg = gasFVF(Pressure, z)
                    Dμ = DgasVis(ρ)
                    Dρ = DgasDensity(Pressure)
                    DBg = DgasFVF(Pressure, z)
                    JE = 1/2*Pre_predict[i,j]*(Dμ*Dρ*Bg + μ*DBg)/(μ*Bg)^2 + (μ*Bg - Pre_predict[i,j+1]/2*(Dμ*Dρ*Bg + μ*DBg))/(μ*Bg)^2
                    J_E = JE * E
                    C_E = (μ*Bg - Pre_predict[i,j+1]/2*(Dμ*Dρ*Bg + μ*DBg))/(μ*Bg)^2*E
                    C_E1 = -1/2*Pre_predict[i,j+1]*(Dμ*Dρ*Bg + μ*DBg)/(μ*Bg)^2*E
                    J[Order[i,j], Order[i,j+1]] = J_E
                    J[Order[i,j], Order[i,j]] -= C_E
                    J[Order[i,j], Order[i,j]] += C_E1
                end
                #C = -(C_N + C_S + C_W + C_E + Gamma[i,j]/Z_predict[i,j])
                J[Order[i,j], Order[i,j]] -= Gamma[i,j]*(gasZ(Pre_predict[i,j]) - Pre_predict[i,j]*DgasZ(Pre_predict[i,j]))/gasZ(Pre_predict[i,j])^2
            end
        end
    end
    return J
end

J = Jacobian_matrix(Pre_predict, Order, Gamma, cnt)
# # The main loop to solve the problem
# Delta_p = J\Residual_matrix
# Pre_var = zeros(cnt,1)
# for i in 1:Row
#     for j in 1:Col
#         if Order[i,j] > 0
#             Pre_var[Order[i,j], 1] = Pre_predict[i,j]
#         end
#     end
# end
# Pre_predict = Pre_var
# Pre_predict = Pre_predict - Delta_p
# println(maximum(Pre_predict))
# println(minimum(Pre_predict))
# count = 0
# for i in 1:length(Pre_predict)
#     if Pre_predict[i] > 0
#         global count += 1
#     end
# end
# println(count)
# println(count == cnt)


iter = 0
while true
    #calculate everything for the loop
    #Initiaze the pressure again
    N,S,W,E,C,Q = construct_Tran_matrix(Pre_init, Pre_predict, Gamma, Order)
    LHS, RHS  = construct_Coefficient_matrix(N,S,W,E,C,Q,Order,cnt)
    Residual_matrix = Residual(LHS, RHS, Pre_predict, Order, cnt)
    J = Jacobian_matrix(Pre_predict, Order, Gamma, cnt)
    Delta_p = J\Residual_matrix
    println(Delta_p[1])
    #Change pressure matrix to vector
    Pre_var = zeros(cnt,1)

    for i in 1:Row
        for j in 1:Col
            if Order[i,j] > 0
                Pre_var[Order[i,j], 1] = Pre_predict[i,j]
            end
        end
    end
    Pre_predict_new = Pre_var
    # for i in 1:Row
    #     for j in 1:Col
    #         if Order[i,j] > 0
    #             Pre_init[i,j] = Pre_var[Order[i,j], 1]
    #         end
    #     end
    # end
    Pre_predict_new = Pre_predict_new - Delta_p
    #println(maximum(Delta_p))
    #println(Pre_predict)
    #change pressure vector to pressure matrix
    Pre_update = zeros(Row,Col)
    for i in 1:Row
        for j in 1:Col
            if Order[i,j] > 0
                Pre_update[i,j] = Pre_predict_new[Order[i,j], 1]
            end
        end
    end
    if maximum(broadcast(abs, Pre_predict_new-Pre_var)) < 0.01
        break
    end

    global Pre_predict = Pre_update
    #Update the pressure
    #Pre_init, Pre_predict = Pre_predict, Pre_update
    #Update Z
    #Z_predict_local = gasCompressibility(Pre_update, 150., Order)
    #global Z_predict = Z_predict_local
    global iter += 1
    #println(maximum(Pre_predict))
    #println(Pre_predict)

end
# println(iter)
# println((Pre_predict[9,11]))
# println(J[50,:])
