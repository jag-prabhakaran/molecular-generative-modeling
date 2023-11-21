

const VisualOutput = ({ VisualData, SMILES_LIST }) => {
    return (
        <div>
            {VisualData ? (
            <div class="VisOutputText" style={{
                fontFamily: "monospace",
                position: "absolute",
                top: "60vh",
                left: "30vw",
                fontSize: "30px",
                color: "black"
            }}>
                Generated Molecules
            </div>
            ): null}
            <diV style= {{
                top: "70vh",
                left: "30vw"
            }}>
                SMILES_LIST
            </diV>
        </div>

    );
};

export default VisualOutput;