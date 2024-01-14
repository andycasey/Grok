
using Korg

function precompile()
    atm = Korg.interpolate_marcs(5777, 4.44)
    ll = Korg.get_APOGEE_DR17_linelist()
    A_X = Korg.format_A_X(0)
    synth = Korg.synthesize(atm, ll, A_X, 15000, 15010, 0.01)
    true
end

#precompile()

function main()
    bytes = 0
    while !eof(stdin)
        line = readline(stdin)
        result = eval(Meta.parse(line))
        println("$line>>>$result")
        bytes += length(line)
    end
    println(bytes)
end

main()

