import { Jsme } from "jsme-react";

const JSMEContainer = ({ onSmilesChange }) => {
    const handleSmilesChange = (newSmiles) => {
        onSmilesChange(newSmiles)
    };

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
            onChange={handleSmilesChange}
            />
    </div>
    );
};

export default JSMEContainer;
