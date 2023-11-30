import React from "react";

const ModelDescription = ({ description, loading }) => {
    return (
        <div>
            <div className = "model-description" style = {{
                position: 'absolute',
                width: '400px',
                height: '525px',
                top: '19vh',
                left: '2vw',
                backgroundColor: '#D0C1C1',
                borderRadius: '5px'
            }}>
                <div className = "model-des-text" style = {{
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
                    <div>
         {loading && <svg style={{
            width:"200px",
            height:"200px",}}>
<defs>
<clipPath id="ldio-ml4ppezh09m-cp" clipPathUnits="userSpaceOnUse">
<rect x="0" y="50" width="100" height="50"></rect>
</clipPath>
<pattern id="ldio-ml4ppezh09m-pattern" patternUnits="userSpaceOnUse" x="0" y="0" width="100" height="100">
<rect x="0" y="0" width="100" height="100" fill="#53ab8b"></rect><circle cx="84" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 108;0 -8" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.4545454545454546s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="47" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 132;0 -32" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.1515151515151516s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="21" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 128;0 -28" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.9696969696969697s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="64" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 146;0 -46" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.0303030303030303s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="90" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 110;0 -10" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.5454545454545454s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="70" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 132;0 -32" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.8181818181818181s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="71" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 135;0 -35" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.3939393939393939s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="34" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 147;0 -47" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.6363636363636364s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="5" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 133;0 -33" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.272727272727273s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="59" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 116;0 -16" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.9090909090909091s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="32" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 132;0 -32" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.878787878787879s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="27" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 106;0 -6" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.5454545454545454s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="1" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 122;0 -22" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.727272727272727s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="51" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 125;0 -25" keyTimes="0;1" dur="3.0303030303030303s" begin="-3.0303030303030303s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="8" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 123;0 -23" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.3333333333333333s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="77" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 124;0 -24" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.15151515151515152s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="12" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 129;0 -29" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.1818181818181819s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="78" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 104;0 -4" keyTimes="0;1" dur="3.0303030303030303s" begin="-3s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="41" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 112;0 -12" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.303030303030303s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="9" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 114;0 -14" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.8181818181818183s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="19" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 138;0 -38" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.8181818181818182s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="48" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 116;0 -16" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.7878787878787878s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="9" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 136;0 -36" keyTimes="0;1" dur="3.0303030303030303s" begin="-2s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="94" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 128;0 -28" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.0606060606060606s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="72" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 145;0 -45" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.8181818181818183s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="11" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 147;0 -47" keyTimes="0;1" dur="3.0303030303030303s" begin="-3s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="29" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 113;0 -13" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.272727272727273s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="52" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 132;0 -32" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.7878787878787878s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="76" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 142;0 -42" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.606060606060606s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="97" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 121;0 -21" keyTimes="0;1" dur="3.0303030303030303s" begin="-2s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="45" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 119;0 -19" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.9696969696969697s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="97" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 112;0 -12" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.45454545454545453s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="22" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 135;0 -35" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.5757575757575757s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="51" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 153;0 -53" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.757575757575758s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="9" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 133;0 -33" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.0606060606060606s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="85" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 122;0 -22" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.6666666666666667s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="74" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 147;0 -47" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.09090909090909091s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="30" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 112;0 -12" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.5757575757575758s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="87" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 130;0 -30" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.0606060606060606s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="26" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 109;0 -9" keyTimes="0;1" dur="3.0303030303030303s" begin="-0.12121212121212122s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="69" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 141;0 -41" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.272727272727273s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="25" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 107;0 -7" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.6363636363636365s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="34" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 124;0 -24" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.5757575757575757s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="63" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 132;0 -32" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.757575757575758s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="39" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 128;0 -28" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.8181818181818181s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="59" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 142;0 -42" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.1818181818181819s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="13" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 138;0 -38" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.2121212121212122s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="72" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 127;0 -27" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.6363636363636362s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="91" cy="0" r="2" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 113;0 -13" keyTimes="0;1" dur="3.0303030303030303s" begin="-2.606060606060606s" repeatCount="indefinite"></animateTransform>
</circle><circle cx="70" cy="0" r="3" fill="#82dbb8">
<animateTransform attributeName="transform" type="translate" values="0 139;0 -39" keyTimes="0;1" dur="3.0303030303030303s" begin="-1.7878787878787878s" repeatCount="indefinite"></animateTransform>
</circle></pattern></defs>
<path fill="url(#ldio-ml4ppezh09m-pattern)" clip-path="url(#ldio-ml4ppezh09m-cp)" d="M59,37.3V18.9c3-0.8,5.1-3.6,5.1-6.8C64.1,8.2,61,5,57.1,5H42.9c-3.9,0-7.1,3.2-7.1,7.1c0,3.2,2.2,6,5.1,6.8v18.4c-11.9,3.8-20.6,15-20.6,28.2C20.4,81.8,33.7,95,50,95s29.6-13.2,29.6-29.6C79.6,52.2,70.9,41.1,59,37.3z"></path>
<path fill="none" stroke="#479f84" stroke-width="5" d="M59,37.3V18.9c3-0.8,5.1-3.6,5.1-6.8C64.1,8.2,61,5,57.1,5H42.9c-3.9,0-7.1,3.2-7.1,7.1c0,3.2,2.2,6,5.1,6.8v18.4c-11.9,3.8-20.6,15-20.6,28.2C20.4,81.8,33.7,95,50,95s29.6-13.2,29.6-29.6C79.6,52.2,70.9,41.1,59,37.3z"></path>

</svg>}
          </div>
                </div>
            </div>
        </div>
    );
};

export default ModelDescription