from email.mime import text
import os
import traceback
from flask import (
    Flask,
    render_template,
    render_template_string,
    request,
    redirect,
    url_for,
    session,
    jsonify,
    flash,
    Response,
)
from dotenv import load_dotenv
from supabase import create_client
from supabase_auth import User
from image import smiles_to_svg
import pandas as pd
import google.genai as genai
import match_compounds as match_compounds

# FILTERED_PATH = "datasets/compounds_filtered.csv"
# df = pd.read_csv(FILTERED_PATH)
# df = df[["Name", "Smiles", "Molecular Weight", "AlogP", "HBA", "HBD"]]
# compounds = df.to_dict(orient="records")

load_dotenv()
# Configuring Gemini API
client = genai.Client(api_key=os.environ["GEMINI_API_KEY"])

# load keys
SUPABASE_URL = os.environ["SUPABASE_URL"]
SUPABASE_SERVICE_ROLE_KEY = os.environ["SUPABASE_SERVICE_ROLE_KEY"]
SUPABASE_ANON_KEY = os.environ["SUPABASE_ANON_KEY"]
FLASK_SECRET_KEY = os.environ["FLASK_SECRET_KEY"]

# flask session
app = Flask(__name__)
app.secret_key = FLASK_SECRET_KEY


# helpers for supabase
def sb_public():
    return create_client(SUPABASE_URL, SUPABASE_ANON_KEY)


def sb_service():
    return create_client(SUPABASE_URL, SUPABASE_SERVICE_ROLE_KEY)


def sb_user():
    sb = create_client(SUPABASE_URL, SUPABASE_ANON_KEY)
    sb.auth.set_session(session["access_token"], session["refresh_token"])
    return sb


def is_logged_in():
    return "access_token" in session and "refresh_token" in session


def require_login():
    if not is_logged_in():
        print(f"DEBUG: require_login failed. Session keys: {list(session.keys())}")
        return redirect(url_for("login"))
    return None


def get_user_id():
    u = sb_user().auth.get_user()
    return u.user.id


def save_match(
    sb_client,
    user_id: str,
    compound_name: str,
    smiles: str,
    properties: dict | None = None,
):
    data = {
        "user_id": user_id,
        "compound_name": compound_name,
        "smiles": smiles,
    }

    if properties:
        data.update(
            {
                "molecular_weight": str(properties.get("molecular_weight", "")),
                "lipophilicity": str(properties.get("lipophilicity", "")),
                "hydrogen_bonding_acceptors": str(
                    properties.get("hydrogen_bonding_acceptors", "")
                ),
                "hydrogen_bonding_donors": str(
                    properties.get("hydrogen_bonding_donors", "")
                ),
            }
        )

    check_res = (
        sb_service()
        .table("matches")
        .select("compound_name")
        .eq("compound_name", data["compound_name"])
        .execute()
    )
    compound_exists = len(check_res.data) > 0

    if not compound_exists:
        sb_client.table("matches").insert(data).execute()


@app.route("/", methods=["GET"])
def root():
    return render_template("landing.html")


@app.route("/onboarding", methods=["GET", "POST"])
def onboarding():
    guard = require_login()
    if guard:
        return guard
    if request.method == "POST":
        try:
            molw = request.form.get("molw")
            lip = request.form.get("lip")
            hba = request.form.get("hba")
            hbd = request.form.get("hbd")

            user_id = get_user_id()

            update_data = {
                "user_id": user_id,
                "molecular_weight": molw,
                "lipophilicity": lip,
                "hydrogen_bonding_acceptors": hba,
                "hydrogen_bonding_donors": hbd,
                "index": 0,
            }

            # Check if user exists in prefs
            try:
                check_res = (
                    sb_service()
                    .table("prefs")
                    .select("user_id")
                    .eq("user_id", user_id)
                    .execute()
                )
                user_exists = len(check_res.data) > 0
            except Exception as e:
                user_exists = False

            if user_exists:
                response = (
                    sb_service()
                    .table("prefs")
                    .update(update_data)
                    .eq("user_id", user_id)
                    .execute()
                )
            else:
                response = sb_service().table("prefs").insert(update_data).execute()

        except Exception as e:
            flash("Something went wrong. Please try again.", "error")
            return render_template("onboarding.html")
        return redirect(url_for("matching"))
    return render_template("onboarding.html")


@app.route("/matching", methods=["GET"])
def matching():
    guard = require_login()
    if guard:
        return guard

    prefs = get_preferences()

    compounds = match_compounds(
        prefs["molecular_weight"],
        prefs["lipophilicity"],
        prefs["hydrogen_bonding_acceptors"],
        prefs["hydrogen_bonding_donors"],
    )

    if "index" not in session:
        session["index"] = 0
        session["accepted"] = []

    idx = prefs["index"]

    if idx >= len(compounds):
        return render_template("matching.html", done=True)

    compound = compounds[idx]

    return render_template(
        "matching.html",
        name=compound.get("Name"),
        smiles=compound.get("Smiles"),
        molecular_weight=compound.get("Molecular Weight"),
        lipophilicity=compound.get("AlogP"),
        hba=compound.get("HBA"),
        hbd=compound.get("HBD"),
        done=False,
    )


@app.route("/api/match", methods=["POST"])
def api_match():
    guard = require_login()
    if guard:
        return jsonify({"error": "unauthorized"}), 401

    payload = request.get_json()
    action = payload.get("action")

    prefs = get_preferences()

    compounds = match_compounds(
        prefs["molecular_weight"],
        prefs["lipophilicity"],
        prefs["hydrogen_bonding_acceptors"],
        prefs["hydrogen_bonding_donors"],
    )

    idx = prefs["index"]

    if idx >= len(compounds):
        return jsonify({"done": True})

    compound = compounds[idx]
    user_id = get_user_id()

    if action == "like":
        properties = {
            "molecular_weight": compound.get("Molecular Weight"),
            "lipophilicity": compound.get("AlogP"),
            "hydrogen_bonding_acceptors": compound.get("HBA"),
            "hydrogen_bonding_donors": compound.get("HBD"),
        }

        print("DEBUG: Saving match:", compound["Name"])

        save_match(
            sb_user(),
            user_id,
            compound.get("Name"),
            compound.get("Smiles"),
            properties,
        )

    update_index(idx + 1)

    next_compound = compounds[idx + 1]
    return jsonify(
        {
            "done": False,
            "name": next_compound["Name"],
            "smiles": next_compound["Smiles"],
        }
    )


# def get_index():
#     res = (
#         sb_service()
#         .table("prefs")
#         .select("index")
#         .eq("user_id", get_user_id())
#         .single()
#         .execute()
#     )
#     return res.data["index"]


def update_index(new_index):
    return (
        sb_service()
        .table("prefs")
        .update({"index": new_index})
        .eq("user_id", get_user_id())
        .execute()
    )


@app.route("/compound-image", methods=["GET"], endpoint="get_compound_image")
def get_compound_image():
    smiles = request.args.get("smiles", "")
    try:
        return Response(smiles_to_svg(smiles), mimetype="image/svg+xml")
    except Exception:
        return "Invalid SMILES", 400


@app.route("/profile", methods=["GET"])
def profile():
    guard = require_login()
    if guard:
        return guard
    if request.method == "GET":
        try:
            prefResponse = (
                sb_service()
                .table("prefs")
                .select("*")
                .eq("user_id", get_user_id())
                .execute()
            )

            print(prefResponse)

            matchResponse = (
                sb_service()
                .table("matches")
                .select("*")
                .eq("user_id", get_user_id())
                .execute()
            )

            return render_template(
                "profile.html", prefs=prefResponse.data, matches=matchResponse.data
            )
        except Exception as e:
            print(e)
            return render_template("profile.html")
    return render_template("profile.html")


def get_preferences():
    response = (
        sb_service().table("prefs").select("*").eq("user_id", get_user_id()).execute()
    )

    res = {
        "molecular_weight": "n/a",
        "lip": "n/a",
        "hba": "n/a",
        "hbd": "n/a",
        "index": 0,
    }

    if response.data:
        res = response.data[0]

    return res


@app.route("/api/profile-chat", methods=["POST"])
def profile_chat():
    guard = require_login()
    if guard:
        return guard

    prompt_type = request.json.get("action")
    user_id = get_user_id()

    print(f"DEBUG: Prompt received: {prompt_type}, user_id: {user_id}")

    prefs_res = sb_service().table("prefs").select("*").eq("user_id", user_id).execute()
    matches_res = (
        sb_service().table("matches").select("*").eq("user_id", user_id).execute()
    )

    print(f"DEBUG: prefs_res: {prefs_res.data}, matches_res: {matches_res.data}")
    prefs = prefs_res.data[0] if prefs_res.data else {}
    matches = matches_res.data

    if not matches:
        return jsonify({"reply": "You have not matched with any compounds yet!"})

    context = {
        "preferences": prefs,
        "matches": [
            {
                "name": m["compound_name"],
                "molecular_weight": m["molecular_weight"],
                "lipophilicity": m["lipophilicity"],
                "hba": m["hydrogen_bonding_acceptors"],
                "hbd": m["hydrogen_bonding_donors"],
            }
            for m in matches
        ],
    }

    system_prompt = f"""
            You are a playful scientific personality analyst.
            You ONLY use the data provided.
            Do NOT invent compounds or preferences.

        User data:
        {context}

        User request:
        {prompt_type}

        Respond in a fun but concise card-style explanation, 
        like what makes the user unique scientifically, 
        based on their compound matches and preferences, 
    no more than 150 words.
"""

    response = client.models.generate_content(
        model="gemini-2.5-flash", contents=system_prompt
    )
    reply_text = response.text
    print(f"DEBUG: Gemini response: {response}")
    return jsonify({"reply": reply_text})


# auth routes
@app.get("/signup")
def signup():
    return render_template("signup.html")


@app.post("/signup")
def signup_post():  # create user in supabase
    email = request.form.get("email", "").strip()
    password = request.form.get("password", "")

    if not email or not password:
        flash("Email and password required.")
        return redirect(url_for("signup"))

    try:
        resp = sb_public().auth.sign_up({"email": email, "password": password})
    except Exception as e:
        flash(f"Signup failed: {e}")
        return redirect(url_for("signup"))  # error

    if getattr(resp, "session", None) is None:
        flash("Account created. Now log in (and confirm email if required).")
        return redirect(url_for("login"))  # for authentictaed login later

    # store token in flask session cookies
    session["access_token"] = resp.session.access_token
    session["refresh_token"] = resp.session.refresh_token
    # redirect to onboarding
    return redirect(url_for("onboarding"))


@app.get("/login")
def login():
    return render_template("login.html")


@app.post("/login")
def login_post():
    email = request.form.get("email", "").strip()
    password = request.form.get("password", "")

    # authenticate
    if not email or not password:
        flash("Email and password required.")
        return redirect(url_for("login"))
    try:
        resp = sb_public().auth.sign_in_with_password(
            {"email": email, "password": password}
        )
    except Exception:
        flash("Login failed. Check email/password.")
        return redirect(url_for("login"))

    session["access_token"] = resp.session.access_token
    session["refresh_token"] = resp.session.refresh_token
    return redirect(url_for("profile"))


@app.post("/logout")
def logout():
    session.clear()
    return redirect(url_for("login"))
