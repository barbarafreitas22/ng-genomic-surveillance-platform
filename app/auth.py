from __future__ import annotations

import os
from pathlib import Path

import streamlit as st
import streamlit_authenticator as stauth
import yaml
from yaml.loader import SafeLoader

_USERS_FILE = Path(__file__).parent / "users.yaml"
_COOKIE_NAME = "ng_surveillance_auth"
_COOKIE_KEY = os.environ.get("AUTH_COOKIE_KEY", "kY16g&wLsIzfeS26#arx^2%*#OAcAldl")
_COOKIE_EXPIRY_DAYS = 30


def _load_credentials() -> dict:
    if not _USERS_FILE.exists():
        return {"usernames": {}}
    with open(_USERS_FILE) as f:
        data = yaml.load(f, Loader=SafeLoader)
    return data.get("credentials", {"usernames": {}})


def get_authenticator() -> stauth.Authenticate:
    return stauth.Authenticate(
        _load_credentials(),
        _COOKIE_NAME,
        _COOKIE_KEY,
        _COOKIE_EXPIRY_DAYS,
        auto_hash=False,
    )


def require_login(authenticator: stauth.Authenticate) -> None:
    if not st.session_state.get("authentication_status"):
        authenticator.login(location="main")
        if st.session_state.get("authentication_status") is False:
            st.error("Username or password incorrect.")
        st.stop()
