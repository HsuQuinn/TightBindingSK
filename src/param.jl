module Params
using LinearAlgebra
#export get_hop_int
#[s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2, S]
"""
Get hopping integrals for tight-binding models.
"""
function hsoc(o1,o2,s1,s2)
    if o1 == "px" && s1 == "up" && o2 == "py" && s2 == "up"
        return -1.0im
    elseif o1 == "px" && s1 == "up" && o2 == "pz" && s2 == "down"
        return 1.0
    elseif o1 == "px" && s1 == "down" && o2 == "py" && s2 == "down"
        return 1.0im 
    elseif o1 == "px" && s1 == "down" && o2 == "pz" && s2 == "up"
        return -1.0
    elseif o1 == "py" && s1 == "up" && o2 == "px" && s2 == "up"
        return 1.0im 
    elseif o1 == "py" && s1 == "up" && o2 == "pz" && s2 == "down"
        return -1.0im 
    elseif o1 == "py" && s1 == "down" && o2 == "px" && s2 == "down"
        return -1.0im 
    elseif o1 == "py" && s1 == "down" && o2 == "pz" && s2 == "up"
        return -1.0im
    elseif o1 == "pz" && s1 == "up" && o2 == "px" && s2 == "down"
        return -1.0
    elseif o1 == "pz" && s1 == "up" && o2 == "py" && s2 == "down"
        return 1.0im 
    elseif o1 == "pz" && s1 == "down" && o2 == "px" && s2 == "up"
        return 1.0
    elseif o1 == "pz" && s1 == "down" && o2 == "py" && s2 == "up"
        return 1.0im 
    end
    return 0
end

function onsite(orbtype, interactions)
    e_s = get(interactions["i"], "e_s", 0.0)
    e_p = get(interactions["i"], "e_p", 0.0)
    e_d = get(interactions["i"], "e_d", 0.0)
    e_S = get(interactions["i"], "e_S", 0.0)
    first_char = first(orbtype)
    if first_char == 's'
        return e_s
    elseif first_char == 'p'
        return e_p
    elseif first_char == 'd'
        return e_d
    elseif first_char == 'S'
        return e_S 
    else
        error("Unsupported orbital type: $orbtype")
    end
end
function get_hop_int(interactions::Dict, l::Float64, m::Float64, n::Float64)
    V_sss = get(interactions["ii"], "V_sss", 0.0)
    V_sps = get(interactions["ii"], "V_sps", 0.0)
    V_pps = get(interactions["ii"], "V_pps", 0.0)
    V_ppp = get(interactions["ii"], "V_ppp", 0.0)
    V_sds = get(interactions["ii"], "V_sds", 0.0)
    V_pds = get(interactions["ii"], "V_pds", 0.0)
    V_pdp = get(interactions["ii"], "V_pdp", 0.0)
    V_dds = get(interactions["ii"], "V_dds", 0.0)
    V_ddp = get(interactions["ii"], "V_ddp", 0.0)
    V_ddd = get(interactions["ii"], "V_ddd", 0.0)
    V_SSs = get(interactions["ii"], "V_SSs", 0.0)
    V_sSs = get(interactions["ii"], "V_sSs", 0.0)
    V_Sps = get(interactions["ii"], "V_Sps", 0.0)
    V_Sds = get(interactions["ii"], "V_Sds", 0.0)

    hop_int = zeros(10,10)

    hop_int[1, 1] = V_sss
    hop_int[1, 2] = l * V_sps
    hop_int[1, 3] = m * V_sps
    hop_int[1, 4] = n * V_sps
    hop_int[2, 1] = -hop_int[1, 2]
    hop_int[3, 1] = -hop_int[1, 3]
    hop_int[4, 1] = -hop_int[1, 4]
    hop_int[1, 5] = sqrt(3) * l * m * V_sds
    hop_int[1, 6] = sqrt(3) * m * n * V_sds
    hop_int[1, 7] = sqrt(3) * l * n * V_sds
    hop_int[5, 1] = hop_int[1, 5]
    hop_int[6, 1] = hop_int[1, 6]
    hop_int[7, 1] = hop_int[1, 7]
    hop_int[1, 8] = sqrt(3) / 2 * (l^2 - m^2) * V_sds
    hop_int[8, 1] = hop_int[1, 8]
    hop_int[1, 9] = (n^2 - 0.5 * (l^2 + m^2)) * V_sds
    hop_int[9, 1] = hop_int[1, 9]

    hop_int[2, 2] = l^2 * V_pps + (1 - l^2) * V_ppp
    hop_int[2, 3] = l * m * (V_pps - V_ppp)
    hop_int[3, 2] = hop_int[2, 3]
    hop_int[2, 4] = l * n * (V_pps - V_ppp)
    hop_int[4, 2] = hop_int[2, 4]

    hop_int[2, 5] = sqrt(3) * l^2 * m * V_pds + m * (1 - 2 * l^2) * V_pdp
    hop_int[2, 6] = l * m * n * (sqrt(3) * V_pds - 2 * V_pdp)
    hop_int[2, 7] = sqrt(3) * l^2 * n * V_pds + n * (1 - 2 * l^2) * V_pdp
    hop_int[5, 2] = -hop_int[2, 5]
    hop_int[6, 2] = -hop_int[2, 6]
    hop_int[7, 2] = -hop_int[2, 7]

    hop_int[2, 8] = 0.5 * sqrt(3) * l * (l^2 - m^2) * V_pds + l * (1 - l^2 + m^2) * V_pdp
    hop_int[2, 9] = l * (n^2 - 0.5 * (l^2 + m^2)) * V_pds - sqrt(3) * l * n^2 * V_pdp
    hop_int[8, 2] = -hop_int[2, 8]
    hop_int[9, 2] = -hop_int[2, 9]

    hop_int[3, 3] = m^2 * V_pps + (1 - m^2) * V_ppp
    hop_int[3, 4] = m * n * (V_pps - V_ppp)
    hop_int[4, 3] = hop_int[3, 4]

    hop_int[3, 5] = sqrt(3) * m^2 * l * V_pds + l * (1 - 2 * m^2) * V_pdp
    hop_int[3, 6] = sqrt(3) * m^2 * n * V_pds + n * (1 - 2 * m^2) * V_pdp
    hop_int[3, 7] = l * m * n * (sqrt(3) * V_pds - 2 * V_pdp)
    hop_int[5, 3] = -hop_int[3, 5]
    hop_int[6, 3] = -hop_int[3, 6]
    hop_int[7, 3] = -hop_int[3, 7]

    hop_int[3, 8] = 0.5 * sqrt(3) * m * (l^2 - m^2) * V_pds - m * (1 + l^2 - m^2) * V_pdp
    hop_int[3, 9] = m * (n^2 - 0.5 * (l^2 + m^2)) * V_pds - sqrt(3) * m * n^2 * V_pdp
    hop_int[8, 3] = -hop_int[3, 8]
    hop_int[9, 3] = -hop_int[3, 9]

    hop_int[4, 4] = n^2 * V_pps + (1 - n^2) * V_ppp

    hop_int[4,5] = l * m * n * (sqrt(3) * V_pds - 2 * V_pdp)
    hop_int[4,6] = sqrt(3) * n^2 * m * V_pds + m * (1.0 - 2 * n^2) * V_pdp
    hop_int[4,7] = sqrt(3) * n^2 * l * V_pds + l * (1.0 - 2 * n^2) * V_pdp
    hop_int[5,4] = -hop_int[4,5]
    hop_int[6,4] = -hop_int[4,6]
    hop_int[7,4] = -hop_int[4,7]

    hop_int[4,8] = 0.5 * sqrt(3) * n * (l^2 - m^2) * V_pds - n * (l^2 - m^2) * V_pdp
    hop_int[4,9] = n * (n^2 - 0.5*(l^2 + m^2)) * V_pds + sqrt(3) * n * (l^2 + m^2) * V_pdp
    hop_int[8,4] = -hop_int[4,8]
    hop_int[9,4] = -hop_int[4,9]

    hop_int[5,5] = l^2 * m^2 * (3 * V_dds - 4 * V_ddp + V_ddd) + (l^2 + m^2) * V_ddp + n^2 * V_ddd
    hop_int[6,6] = m^2 * n^2 * (3 * V_dds - 4 * V_ddp + V_ddd) + (m^2 + n^2) * V_ddp + l^2 * V_ddd
    hop_int[7,7] = n^2 * l^2 * (3 * V_dds - 4 * V_ddp + V_ddd) + (n^2 + l^2) * V_ddp + m^2 * V_ddd

    hop_int[5,6] = l * m^2 * n * (3 * V_dds - 4 * V_ddp + V_ddd) + l * n * (V_ddp - V_ddd)
    hop_int[5,7] = n * l^2 * m * (3 * V_dds - 4 * V_ddp + V_ddd) + n * m * (V_ddp - V_ddd)
    hop_int[6,7] = m * n^2 * l * (3 * V_dds - 4 * V_ddp + V_ddd) + m * l * (V_ddp - V_ddd)
    hop_int[6,5] = hop_int[5,6]
    hop_int[7,5] = hop_int[5,7]
    hop_int[7,6] = hop_int[6,7]

    hop_int[5,8] = 0.5 * l * m * (l^2 - m^2) * (3 * V_dds - 4 * V_ddp + V_ddd)
    hop_int[6,8] = 0.5 * m * n * ((l^2 - m^2) * (3 * V_dds - 4 * V_ddp + V_ddd) - 2 * (V_ddp - V_ddd))
    hop_int[7,8] = 0.5 * n * l * ((l^2 - m^2) * (3 * V_dds - 4 * V_ddp + V_ddd) + 2 * (V_ddp - V_ddd))
    hop_int[8,5] = hop_int[5,8]
    hop_int[8,6] = hop_int[6,8]
    hop_int[8,7] = hop_int[7,8]

    hop_int[5,9] = sqrt(3) * (l * m * (n^2 - 0.5*(l^2 + m^2)) * V_dds - 2 * l * m * n^2 * V_ddp + 0.5 * l * m * (1.0 + n^2) * V_ddd)
    hop_int[6,9] = sqrt(3) * (m * n * (n^2 - 0.5*(l^2 + m^2)) * V_dds + m * n * (l^2 + m^2 - n^2) * V_ddp - 0.5 * m * n * (l^2 + m^2) * V_ddd)
    hop_int[7,9] = sqrt(3) * (n * l * (n^2 - 0.5*(l^2 + m^2)) * V_dds + n * l * (l^2 + m^2 - n^2) * V_ddp - 0.5 * n * l * (l^2 + m^2) * V_ddd)
    hop_int[9,5] = hop_int[5,9]
    hop_int[9,6] = hop_int[6,9]
    hop_int[9,7] = hop_int[7,9]

    hop_int[8,8] = 0.25 * (l^2 - m^2)^2 * (3 * V_dds - 4 * V_ddp + V_ddd) + (l^2 + m^2) * V_ddp + n^2 * V_ddd
    hop_int[9,9] = 0.75 * (l^2 + m^2)^2 * V_ddd + 3 * (l^2 + m^2) * n^2 * V_ddp + 0.25 * (l^2 + m^2 - 2 * n^2)^2 * V_dds
    hop_int[8,9] = sqrt(3) * 0.25 * (l^2 - m^2) * (n^2 * (2 * V_dds - 4 * V_ddp + V_ddd) + V_ddd - (l^2 + m^2) * V_dds)
    hop_int[9,8] = sqrt(3) * 0.25 * (l^2 - m^2) * (n^2 * (2 * V_dds - 4 * V_ddp + V_ddd) + V_ddd - (l^2 + m^2) * V_dds)

    hop_int[10,10] = V_SSs
    hop_int[1,10] = V_sSs
    hop_int[10,1] = V_sSs
    hop_int[10,2] = l * V_Sps
    hop_int[10,3] = m * V_Sps
    hop_int[10,4] = n * V_Sps
    hop_int[2,10] = -hop_int[10,2]
    hop_int[3,10] = -hop_int[10,3]
    hop_int[4,10] = -hop_int[10,4]
    hop_int[10,5] = sqrt(3) * l * m * V_Sds
    hop_int[10,6] = sqrt(3) * m * n * V_Sds
    hop_int[10,7] = sqrt(3) * l * n * V_Sds
    hop_int[5,10] = hop_int[10,5]
    hop_int[6,10] = hop_int[10,6]
    hop_int[7,10] = hop_int[10,7]
    hop_int[10,8] = (sqrt(3)/2.0) * (l * l - m * m) * V_Sds
    hop_int[8,10] = hop_int[10,8]
    hop_int[10,9] = (n^2 - 0.5*(l^2 + m^2)) * V_Sds
    hop_int[9,10] = hop_int[10,9]
    return hop_int
end

end # module Params