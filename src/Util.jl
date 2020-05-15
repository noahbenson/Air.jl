################################################################################
# Util.jl
# Utilities for Air that don't depend on other components of Air.
# by Noah C. Benson

# This helper function makes sure than a function-arg declaration has a name.
# This helper function makes sure than a function-arg declaration has a name.
_memoize_fixarg(arg::Expr) = (arg.head == :(::) && length(arg.args) == 1
                              ? Expr(:(::), gensym(), arg.args[1])
                              : arg)
"""
    @memoize name(args...) = expr
    @memoize name(args...) where {...} = expr

@memoize is a macro for declaring that the function declaration that follows                                                                                                                               should be memoized in a private dictionary and any pre-calculated value should                                                                                                                             be returned from that dictionary instead of being recalculated.
"""
macro memoize(assgn::Expr)
    (assgn.head == :(=)) || throw(
        ArgumentError("memconst must be given an assignment expression"))
    # Parse the assignment statement.
    lhs  = assgn.args[1]
    expr = assgn.args[2]
    if lhs.head == :call
        fsym = lhs.args[1]
        args = [_memoize_fixarg(a) for a in lhs.args[2:end]]
        lhs = Expr(:call, fsym, args...)
    elseif lhs.head == :where
        fsig = lhs.args[1]
        fsym = fsig.args[1]
        args = [_memoize_fixarg(a) for a in fsig.args[2:end]]
        fsig = Expr(:call, fsym, args...)
        lhs = Expr(:where, fsig, lhs.args[2:end]...)
    else
        throw(ArgumentError("memconst assignment LHS must be a call or where expression"))
    end
    # Make an expression for the tuple of arguments.
    argtup = Expr(:tuple, args...)
    # Symbols we will need in the generated code.
    sDict = gensym()
    sLock = gensym()
    sTmp  = gensym()
    quote
        $sLock = ReentrantLock()
        $sDict = Dict{Tuple, Any}()
        $lhs = begin
            lock($sLock)
            $sTmp = get!($sDict, $argtup) do; $expr end
            unlock($sLock)
            return $sTmp
        end
        function forget(::typeof($fsym))
            global $sDict
            $sTmp = $sDict
            $sDict = Dict{Tuple, Any}()
            return $sTmp
        end
        Base.keys(::typeof($fsym)) = keys($sDict)
        $fsym
    end |> esc
end
