"""
    markgraph_todot(mg) -> String

Return a Graphviz DOT description of marking graph `mg`. Each node is a marking,
labelled by its token vector and styled by `StateType`:

- tangible  — solid circle
- vanishing — dashed circle (left in zero time)
- absorbing — double circle

Each edge is labelled with the firing transition's label and its rate (`exp`
edges), weight (`imm` edges), or `gen` for a general-transition jump. The initial
marking is drawn bold.
"""
function markgraph_todot(mg::MarkingGraph)
    io = IOBuffer()
    println(io, "digraph { layout=dot; overlap=false; splines=true; node [fontsize=10];")
    trlabel(id) = mg.pn.trans[id].label
    markstr(i) = "[" * join(mg.states[i], ",") * "]"

    for i in 1:nstates(mg)
        shape = mg.types[i] == ABSORBING ? "doublecircle" : "circle"
        style = mg.types[i] == VANISHING ? ", style=dashed" : ""
        bold = i == mg.initial_state ? ", penwidth=2" : ""
        println(io, "\"s$i\" [shape=$shape, label=\"$(markstr(i))\"$style$bold];")
    end
    for e in mg.edges
        tag = e.kind === :exp ? "$(round(e.value; digits=4))" :
              e.kind === :imm ? "w=$(round(e.value; digits=4))" : "gen"
        println(io, "\"s$(e.src)\"->\"s$(e.dst)\" [label=\"$(trlabel(e.trid)) ($tag)\"];")
    end
    println(io, "}")
    String(take!(io))
end
