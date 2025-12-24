import groovy.json.JsonOutput

process SAVE_PARAMS {

    label 'process_single'

    tag "Save final parameters."

    input:
    val workflowParams

    output:
    path "params.json"

    script:
    def jsonStr  = JsonOutput.toJson(workflowParams)
    def pretty   = JsonOutput.prettyPrint(jsonStr)
    """
    cat << 'EOF' > params.json
    ${pretty}
    EOF
    """
}
