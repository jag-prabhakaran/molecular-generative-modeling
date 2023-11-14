
const GenerateButton = ({ onGenerate }) => {
    return (
        <div>
            <button id="generate" onClick={onGenerate}
            style={{
                width: "243px",
                height: "52px",
                top: "52.5vh",
                left: "30vw",
                position: "absolute",
                background: "#8BCACC",
                fontFamily: "monospace"
            }}
            >
                GENERATE
            </button>
        </div>
    );
};
export default GenerateButton;