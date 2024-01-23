
const Properties = ({ logPMax, onLogPMaxChange, logPMin, onLogPMinChange, NumOfMolecules, onNumOfMoleculesChange, qedMin, onQedMinChange, qedMax, onQedMaxChange }) => {
    return (
        <div style={{
            position: "absolute",
            left: "75vw"
        }}>
            <div className="PropertiesText" style={{
                top: "-39vh",
                color: "black",
                fontSize: "2.5vh",
                fontFamily: "monospace",
                fontWeight: 400,
                wordWrap: "break-word",
                position: "relative"
            }}>
                Properties
            </div>
            <div className="logP" style={{
                background: "#8BCACC",
                borderRadius: "5px",
                width: "300px",
                height: "50px",
                top: "-37vh",
                position: "relative"
            }}>
                <div className="logPText" style={{
                    color: "black",
                    fontSize: "2vh",
                    fontFamily: "monospace",
                    fontWeight: 200,
                    wordWrap: "break-word",
                    top: "1.5vh",
                    left: "1.5vw",
                    position: "absolute"
                }}>
                    logP
                </div>
                <div className="logPMaxValue">
                    <input style= {{
                        position: "absolute",
                        top: "0.7vh",
                        left: "5vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={logPMax}
                    id="logPMax"
                    onChange={(e) => onLogPMaxChange(e.target.value)}>
                    </input>
                </div>
                <div className="ToText" style={{
                    color: "black",
                    fontSize: "2vh",
                    fontFamily: "monospace",
                    fontWeight: "150",
                    wordWrap: "break-word",
                    top: "1.5vh",
                    left: "9vw",
                    position: "absolute"
                }}>
                    to
                </div>
                <div className="logPMinValue">
                    <input style= {{
                        position: "absolute",
                        top: "0.7vh",
                        left: "11vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={logPMin}
                    id="logPMin"
                    onChange={(e) => onLogPMinChange(e.target.value)}
                    >
                    </input>
                </div>
            </div>
            <div className="NumOfMolecules" style={{
                background: "#8BCACC",
                borderRadius: "5px",
                width: "300px",
                height: "50px",
                top: "-35vh",
                position: "relative"
            }}>
                <div className="NumOfMoleculesText" style={{
                    color: "black",
                    fontSize: "2vh",
                    fontFamily: "monospace",
                    fontWeight: 200,
                    wordWrap: "break-word",
                    top: "1.2vh",
                    left: "-2vw",
                    position: "relative"
                }}>
                    Num Of Molecules
                </div>
                <div className="NumOfMoleculesValue">
                    <input style= {{
                        position: "relative",
                        top: "-2.5vh",
                        left: "4vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={NumOfMolecules}
                    id="NumOfMoleculesValue"
                    onChange={(e) => onNumOfMoleculesChange(e.target.value)}
                    >
                    </input>
                </div>
            </div>
            <div className= "qed" style={{
                background: "#8BCACC",
                borderRadius: "5px",
                width: "300px",
                height: "50px",
                top: "-33vh",
                position: "relative"
            }}>
                <div className = "qedText" style={{
                    color: "black",
                    fontSize: "2vh",
                    fontFamily: "monospace",
                    fontWeight: 200,
                    wordWrap: "break-word",
                    top: "1.5vh",
                    left: "1.5vw",
                    position: "absolute"    
                }}>
                    QED Score
                </div>
                <div className="qedMaxValue">
                    <input style= {{
                        position: "absolute",
                        top: "0.7vh",
                        left: "5vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={qedMax}
                    id="qedMax"
                    onChange={(e) => onQedMaxChange(e.target.value)}>
                    </input>
                </div>
                <div className="ToText" style={{
                    color: "black",
                    fontSize: "2vh",
                    fontFamily: "monospace",
                    fontWeight: "150",
                    wordWrap: "break-word",
                    top: "1.5vh",
                    left: "9vw",
                    position: "absolute"
                }}>
                    to
                </div>
                <div className="qedMinValue">
                    <input style= {{
                        position: "absolute",
                        top: "0.7vh",
                        left: "11vw",
                        width: "50px",
                        height: "30px",
                        backgroundColor: "lightgray"
                    }}
                    type="number"
                    value={qedMin}
                    id="qedMin"
                    onChange={(e) => onQedMinChange(e.target.value)}
                    >
                    </input>
                </div>
            </div>
        </div>
    )
}

export default Properties