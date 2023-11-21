import React from "react";

const ModelDescription = ({ description }) => {
    return (
        <div>
            <div class = "model-description" style = {{
                position: 'absolute',
                width: '400px',
                height: '365px',
                top: '19vh',
                left: '2vw',
                backgroundColor: '#D0C1C1',
                borderRadius: '5px'
            }}>
                <div class = "model-des-text" style = {{
                    height: '50px',
                    width: '375px',
                    left: '0.25vw',
                    position: 'relative'
                }}>
                    <p style = {{
                        fontSize: '15px',
                        color: '#000',
                        fontFamily: "monospace",
                        wordWrap: "break-word"
                    }}
                        dangerouslySetInnerHTML={{ __html: description }}
                    />
                </div>
            </div>
        </div>
    );
};

export default ModelDescription