import { Jsme } from "jsme-react";

const JSMEContainer = () => {
    return (
    <div style={{
        position: "absolute",
        top: "11vh",
        left: "30vw"
        }}>
            <Jsme
            height="360px"
            width="680px"
            options="newlook,star"
            jme="startingStructure"
            />
    </div>
    );
};

export default JSMEContainer;
