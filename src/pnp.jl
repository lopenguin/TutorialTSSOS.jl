## Functions demonstrating the use of TSSOS for perspective-n-point
# note that required imports are in `TutorialTSSOS.jl`


"""
    genPnPproblem(N, σ, camK; r=0.2)

Generate a PnP problem with `N` points.

# Arguments
- `N`: number of points
- `σ`: Gaussian measurement noise (3D)
- `camK`: camera calibration matrix
- `r=0.2`: radius of generated shape

# Returns:
- `y`: homogenized pixel detections corrupted by noise
- `b`: 3D points in canonical frame
- `gt`: tuple of ground truth (R,t)
"""
function genPnPproblem(N, σ, camK; r=0.2)
    # gt 3D points
    b = r*randn(3, N)
    b .-= mean(b,dims=2)

    # gt position
    t = randn(3) + [0; 0; 1.5]

    # gt rotation: make sure pointing towards object
    found_R = false
    tn = normalize(t)
    R = I
    for _ = 1:100
        R = randrotation()
        if R[:,3]'*tn > cos(60*π/180.) # within 60 deg.
            found_R = true
            break
        end
    end
    if !found_R
        error("Problem generation failed to find R! Suggest rerunning.")
    end

    # save
    gt = (R, t)

    # convert to pixel measurements
    ϵ = randn(3, N) .* σ'
    y = camK*(R*b .+ t + ϵ)
    y = reduce(hcat, eachcol(y) ./ y[3,:])

    # save
    return y, b, gt
end